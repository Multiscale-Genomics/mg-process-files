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
import subprocess
import shlex

import numpy as np
import h5py

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT, IN
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_INOUT, FILE_OUT, IN  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task  # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on  # pylint: disable=ungrouped-imports

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------


class bedIndexerTool(Tool):
    """
    Tool for running indexers over a BED file for use in the RESTful API
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("BED File Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def bed_feature_length(self, file_bed):
        """
        BED Feature Length

        Function to calcualte the averagte length of a feature in BED file.

        Parameters
        ----------
        file_bed : str
            Location of teh BED file

        Returns
        -------
        average_feature_length : int
            The average length of the features in a BED file.
        """

        total_feature_count = 0
        total_feature_length = 0

        with open(file_bed, 'r') as f_in:
            for line in f_in:
                line = line.strip()
                sline = line.split("\t")

                start = int(sline[1])
                end = int(sline[2])
                length = end-start

                total_feature_count += 1
                total_feature_length += length

        return total_feature_length / total_feature_count

    @task(returns=bool, file_sorted_bed=FILE_IN, file_chrom=FILE_IN,
          file_bb=FILE_OUT, bed_type=IN, isModifier=False)
    def bed2bigbed(self, file_sorted_bed, file_chrom, file_bb, bed_type=None):  # pylint: disable=no-self-use
        """
        BED to BigBed converter

        This uses the ``bedToBigBed`` program binary provided at
        http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
        to perform the conversion from bed to bigbed.

        Parameters
        ----------
        file_sorted_bed : str
            Location of the sorted BED file
        file_chrom : str
            Location of the chrom.size file
        file_bb : str
            Location of the bigBed file

        Example
        -------
        .. code-block:: python
           :linenos:

           if not self.bed2bigbed(bed_file, chrom_file, bb_file):
               output_metadata.set_exception(
                   Exception(
                       "bed2bigbed: Could not process files {}, {}.".format(*input_files)))

        """
        command_line = 'bedToBigBed'
        if bed_type is not None:
            command_line += ' -type=' + str(bed_type)

        command_line += ' ' + file_sorted_bed + ' ' + file_chrom + ' ' + file_bb + '.tmp.bb'

        logger.info('BED 2 BIGBED:', command_line)

        args = shlex.split(command_line)
        process_handle = subprocess.Popen(args)
        process_handle.wait()

        with open(file_bb, 'wb') as f_out:
            with open(file_bb + '.tmp.bb', 'rb') as f_in:
                f_out.write(f_in.read())

        return True

    @task(returns=bool, file_id=IN, assembly=IN, file_sorted_bed=FILE_IN,
          file_hdf5=FILE_INOUT)
    def bed2hdf5(self, file_id, assembly, file_sorted_bed, file_hdf5):  # pylint: disable=no-self-use, too-many-locals,too-many-statements,too-many-branches
        """
        BED to HDF5 converter

        Loads the BED file into the HDF5 index file that gets used by the REST
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
        feature_length : int
            Defines the level of resolution that the features should be recorded
            at. The 2 options are 1 or 1000. 1 records features at every single
            base whereas 1000 groups features into 1000bp chunks. The single
            base pair option should really only be used when features are less
            than 10bp to
        file_sorted_bed : str
            Location of the sorted BED file
        file_hdf5 : str
            Location of the HDF5 index file

        Example
        -------
        .. code-block:: python
           :linenos:

           if not self.bed2hdf5(file_id, assembly, bed_file, hdf5_file):
               output_metadata.set_exception(
                   Exception(
                       "bed2hdf5: Could not process files {}, {}.".format(*input_files)))

        """
        max_files = 1024
        max_chromosomes = 1024
        max_chromosome_size = 2000000000

        feature_length = self.bed_feature_length(file_sorted_bed)
        storage_level = 1000
        if feature_length < 10:
            storage_level = 1

        hdf5_in = h5py.File(file_hdf5, "a")

        if str(assembly) in hdf5_in:
            grp = hdf5_in[str(assembly)]
            meta = hdf5_in['meta']

            dset1 = grp['data1']
            dset1k = grp['data1k']
            fset = grp['files']
            cset = grp['chromosomes']
            file_idx_1 = [fs for fs in fset[0] if fs != '']
            file_idx_1k = [fs for fs in fset[1] if fs != '']
            if file_id not in file_idx_1 and file_id not in file_idx_1k:
                if storage_level == 1000:
                    file_idx_1k.append(file_id)
                else:
                    file_idx_1.append(file_id)

                # pylint comment: resize is a valid member of the objects
                dset1.resize((dset1.shape[0], dset1.shape[1] + 1, max_chromosome_size))  # pylint: disable=no-member
                dset1k.resize((dset1k.shape[0], dset1k.shape[1] + 1, max_chromosome_size/1000))  # pylint: disable=no-member
            chrom_idx = [c for c in cset if c != '']

        else:
            # Create the initial dataset with minimum values
            grp = hdf5_in.create_group(str(assembly))
            hdf5_in.create_group('meta')

            dtf = h5py.special_dtype(vlen=str)
            dtc = h5py.special_dtype(vlen=str)
            fset = grp.create_dataset('files', (2, max_files), dtype=dtf)
            cset = grp.create_dataset('chromosomes', (max_chromosomes,), dtype=dtc)

            file_idx_1 = []
            file_idx_1k = []
            chrom_idx = []

            logger.info(str(max_chromosome_size), str(max_chromosomes), str(max_files))
            dset1 = grp.create_dataset(
                'data1', (0, 1, max_chromosome_size),
                maxshape=(max_chromosomes, max_files, max_chromosome_size),
                dtype='bool', chunks=True, compression="gzip"
            )
            dset1k = grp.create_dataset(
                'data1k', (0, 1, max_chromosome_size/1000),
                maxshape=(max_chromosomes, max_files, max_chromosome_size/1000),
                dtype='bool', chunks=True, compression="gzip"
            )

            if storage_level == 1000:
                file_idx_1k.append(file_id)
            else:
                file_idx_1.append(file_id)

        # Save the list of files
        fset[0, 0:len(file_idx_1)] = file_idx_1
        fset[1, 0:len(file_idx_1k)] = file_idx_1k

        file_chrom_count = 0

        if storage_level == 1000:
            dnp = np.zeros([max_chromosome_size/1000], dtype='bool')
        else:
            dnp = np.zeros([max_chromosome_size], dtype='bool')

        previous_chrom = ''

        loaded = False

        with open(file_sorted_bed, 'r') as f_in:
            for line in f_in:
                line = line.strip()
                length = line.split("\t")

                chrom = str(length[0])
                start = int(length[1])
                end = int(length[2])

                loaded = False

                if chrom != previous_chrom and previous_chrom != '':
                    file_chrom_count += 1
                    if previous_chrom not in chrom_idx:
                        chrom_idx.append(previous_chrom)
                        cset[0:len(chrom_idx)] = chrom_idx
                        dset1.resize(
                            (
                                dset1.shape[0]+1,
                                dset1.shape[1],
                                max_chromosome_size)
                        )
                        dset1k.resize(
                            (
                                dset1k.shape[0]+1,
                                dset1k.shape[1],
                                max_chromosome_size/1000
                            )
                        )

                    loaded = True

                    if storage_level == 1000:
                        dset1k[chrom_idx.index(previous_chrom), file_idx_1k.index(file_id), :] = dnp
                        dnp = np.zeros([max_chromosome_size/1000], dtype='bool')
                    else:
                        dset1[chrom_idx.index(previous_chrom), file_idx_1.index(file_id), :] = dnp
                        dnp = np.zeros([max_chromosome_size], dtype='bool')

                previous_chrom = chrom
                if storage_level == 1000:
                    dnp[(start/1000):(end/1000)+1] = '1'
                else:
                    dnp[start:end+1] = '1'

            if loaded is False:
                if previous_chrom not in chrom_idx:
                    chrom_idx.append(chrom)
                    cset[0:len(chrom_idx)] = chrom_idx
                    dset1.resize((dset1.shape[0] + 1, dset1.shape[1], max_chromosome_size))
                    dset1k.resize((dset1k.shape[0] + 1, dset1k.shape[1], max_chromosome_size/1000))

                if storage_level == 1000:
                    dset1k[
                        chrom_idx.index(previous_chrom)/1000, file_idx_1k.index(file_id), :
                    ] = dnp
                else:
                    dset1[
                        chrom_idx.index(previous_chrom), file_idx_1.index(file_id), :
                    ] = dnp

        hdf5_in.close()

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Function to run the BED file sorter and indexer so that the files can
        get searched as part of the REST API

        Parameters
        ----------
        input_files : list
            bed_file : str
                Location of the sorted bed file
            chrom_size : str
                Location of chrom.size file
            hdf5_file : str
                Location of the HDF5 index file
        metadata : list
            file_id : str
                file_id used to identify the original bed file
            assembly : str
                Genome assembly accession

        Returns
        -------
        list
            bed_file : str
                Location of the sorted bed file
            bb_file : str
                Location of the BigBed file
            hdf5_file : str
                Location of the HDF5 index file

        Example
        -------
        .. code-block:: python
           :linenos:

           import tool

           # Bed Indexer
           b = tool.bedIndexerTool(self.configuration)
           bi, bm = bd.run(
               [bed_file_id, chrom_file_id, hdf5_file_id], [], {'assembly' : assembly}
           )
        """
        bed_type = None
        if "bed_type" in self.configuration:
            bed_type = self.configuration['bed_type']

        results = self.bed2bigbed(
            input_files["bed"], input_files["chrom_file"], output_files["bb_file"], bed_type)
        results = compss_wait_on(results)

        results = self.bed2hdf5(
            input_files['bed'], input_metadata["bed"].meta_data["assembly"],
            input_files["bed"], input_files["hdf5_file"]
        )
        results = compss_wait_on(results)

        output_generated_files = {
            "bb_file": output_files["bb_file"],
            "hdf5_file": input_metadata["bed"].meta_data["assembly"]
        }
        output_metadata = {
            "bb_file": input_files["bed"],
            "hdf5_file": input_metadata["hdf5_file"]
        }

        return (output_generated_files, output_metadata)

# ------------------------------------------------------------------------------
