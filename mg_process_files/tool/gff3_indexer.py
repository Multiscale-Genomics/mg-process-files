"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from __future__ import print_function

import sys
import numpy as np
import h5py
import pysam

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_INOUT, FILE_OUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_INOUT, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------


class gff3IndexerTool(Tool):
    """
    Tool for running indexers over a WIG file for use in the RESTful API
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("GFF3 File Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, file_sorted_gff3=FILE_IN, file_sorted_gz_gff3=FILE_OUT,
          file_gff3_tbi=FILE_OUT)
    def gff32tabix(self, file_sorted_gff3, file_sorted_gz_gff3, file_gff3_tbi):  # pylint: disable=no-self-use,unused-argument
        """
        GFF3 to Tabix

        Compresses the sorted GFF3 file and then uses Tabix to generate an index
        of the GFF3 file.

        Parameters
        ----------
        file_sorted_gff3 : str
            Location of a sorted GFF3 file
        file_sorted_gz_gff3 : str
            Location of the bgzip compressed GFF3 file
        file_gff3_tbi : str
            Location of the Tabix index file

        Example
        -------
        .. code-block:: python
           :linenos:

           if not self.gff32tabix(self, file_sorted_gff3, gz_file, tbi_file):
               output_metadata.set_exception(
                   Exception(
                       "gff32tabix: Could not process files {}, {}.".format(*input_files)))
        """
        pysam.tabix_compress(file_sorted_gff3, file_sorted_gz_gff3)  # pylint: disable=no-member
        pysam.tabix_index(file_sorted_gz_gff3, preset='gff')  # pylint: disable=no-member
        return True

    @task(returns=bool, file_id=IN, assembly=IN, file_sorted_gff3=FILE_IN, file_hdf5=FILE_INOUT)
    def gff32hdf5(self, file_id, assembly, file_sorted_gff3, file_hdf5):  # pylint: disable=no-self-use,too-many-locals,too-many-statements
        """
        GFF3 to HDF5 converter

        Loads the GFF3 file into the HDF5 index file that gets used by the REST
        API to determine if there are files that have data in a given region.
        Overlapping regions are condensed into a single feature block rather
        than maintaining all of the detail of the original bed file.

        Parameters
        ----------
        file_id : str
            The file_id as stored by the DM-API so that it can be used for file
            retrieval later
        assembly : str
            Assembly of the genome that is getting indexed so that the
            chromosomes match
        file_sorted_gff3 : str
            Location of the sorted GFF3 file
        file_hdf5 : str
            Location of the HDF5 index file

        Example
        -------
        .. code-block:: python
           :linenos:

           if not self.gff32hdf5(file_id, assembly, bed_file, hdf5_file):
               output_metadata.set_exception(
                   Exception(
                       "gff32hdf5: Could not process files {}, {}.".format(*input_files)))

        """
        max_files = 1024
        max_chromosomes = 1024
        max_chromosome_size = 2000000000

        f_h5_in = h5py.File(file_hdf5, "a")

        if str(assembly) in f_h5_in:
            grp = f_h5_in[str(assembly)]

            dset = grp['data']
            fset = grp['files']
            cset = grp['chromosomes']
            file_idx = [i for i in fset if i != '']
            if file_id not in file_idx:
                file_idx.append(file_id)
                dset.resize((dset.shape[0], dset.shape[1] + 1, max_chromosome_size))  # pylint: disable=no-member
            chrom_idx = [c for c in cset if c != '']

        else:
            # Create the initial dataset with minimum values
            grp = f_h5_in.create_group(str(assembly))
            f_h5_in.create_group('meta')

            dtf = h5py.special_dtype(vlen=str)
            dtc = h5py.special_dtype(vlen=str)
            fset = grp.create_dataset('files', (max_files,), dtype=dtf)
            cset = grp.create_dataset('chromosomes', (max_chromosomes,), dtype=dtc)

            file_idx = [file_id]
            chrom_idx = []

            dset = grp.create_dataset(
                'data', (0, 1, max_chromosome_size),
                maxshape=(max_chromosomes, max_files, max_chromosome_size),
                dtype='bool', chunks=True, compression="gzip"
            )

        # Save the list of files
        fset[0:len(file_idx)] = file_idx

        file_chrom_count = 0

        dnp = np.zeros([max_chromosome_size], dtype='bool')

        previous_chrom = ''
        loaded = False

        with open(file_sorted_gff3, 'r') as f_in:
            for line in f_in:
                if line[0] == '#':
                    continue

                line = line.strip()
                sline = line.split("\t")

                chrom = str(sline[0])
                start = int(sline[3])
                end = int(sline[4])

                loaded = False

                if chrom != previous_chrom and previous_chrom != '':
                    file_chrom_count += 1
                    if previous_chrom not in chrom_idx:
                        chrom_idx.append(previous_chrom)
                        cset[0:len(chrom_idx)] = chrom_idx
                        dset.resize((dset.shape[0] + 1, dset.shape[1], max_chromosome_size))

                    dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp
                    loaded = True

                    if file_chrom_count == 5:
                        break

                    dnp = np.zeros([max_chromosome_size], dtype='bool')

                previous_chrom = chrom
                dnp[start:end + 1] = 1

            if loaded is False:
                if previous_chrom not in chrom_idx:
                    chrom_idx.append(chrom)
                    cset[0:len(chrom_idx)] = chrom_idx
                    dset.resize((dset.shape[0] + 1, dset.shape[1], max_chromosome_size))

                dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp

        f_h5_in.close()

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Function to run the BED file sorter and indexer so that the files can
        get searched as part of the REST API

        Parameters
        ----------
        input_files : list
            gff3_file : str
                Location of the bed file
            hdf5_file : str
                Location of the HDF5 index file
        meta_data : list
            file_id : str
                file_id used to identify the original bed file
            assembly : str
                Genome assembly accession

        Returns
        -------
        list
            gz_file : str
                Location of the sorted gzipped GFF3 file
            tbi_file : str
                Location of the Tabix index file
            hdf5_file : str
                Location of the HDF5 index file
        """
        # handle error
        results_1 = self.gff32tabix(
            input_files["gff3"], output_files["gz_file"],
            output_files["tbi_file"])
        results_1 = compss_wait_on(results_1)

        results_2 = self.gff32hdf5(
            input_files["gff3"],
            input_metadata["gff3"].meta_data["assembly"],
            input_files["gff3"],
            input_files["hdf5_file"])
        results_2 = compss_wait_on(results_2)

        output_generated_files = {
            "gz_file": output_files["gz_file"],
            "tbi_file": output_files["tbi_file"],
            "hdf5_file": input_metadata["gff3"].meta_data["assembly"]
        }
        output_metadata = {
            "gz_file": input_metadata["gff3"],
            "tbi_file": input_metadata["gff3"],
            "hdf5_file": input_metadata["hdf5_file"]
        }

        return (output_generated_files, output_metadata)

# ------------------------------------------------------------------------------
