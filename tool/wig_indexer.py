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


class wigIndexerTool(Tool):
    """
    Tool for running indexers over a WIG file for use in the RESTful API
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        print("WIG File Indexer")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, file_wig=FILE_IN, file_chrom=FILE_IN, file_bw=FILE_OUT,
          isModifier=False)
    def wig2bigwig(self, file_wig, file_chrom, file_bw):  # pylint: disable=no-self-use
        """
        WIG to BigWig converter

        This uses the ``wigToBigWig`` program binary provided at
        http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
        to perform the conversion from WIG to BigWig.

        Parameters
        ----------
        file_wig : str
            Location of the wig file
        file_chrom : str
            Location of the chrom.size file
        file_bw : str
            Location of the bigWig file

        Example
        -------
        .. code-block:: python
           :linenos:

           if not self.wig2bigwig(wig_file, chrom_file, bw_file):
               output_metadata.set_exception(
                   Exception(
                       "wig2bigWig: Could not process files {}, {}.".format(*input_files)))

        """
        command_line = 'wigToBigWig ' + file_wig + ' ' + file_chrom + ' ' + file_bw + '.tmp.bw'
        args = shlex.split(command_line)
        process = subprocess.Popen(args)
        process.wait()

        logger.info('BIGWIG - COMMAND:', command_line)
        logger.info('BIGWIG - FILES:', file_wig, file_chrom, file_bw)

        with open(file_bw, 'wb') as f_out:
            with open(file_bw + '.tmp.bw', 'rb') as f_in:
                f_out.write(f_in.read())

        return True

    @task(returns=bool, file_id=IN, assembly=IN, file_wig=FILE_IN, file_hdf5=FILE_INOUT)
    def wig2hdf5(self, file_id, assembly, file_wig, file_hdf5):  # pylint: disable=no-self-use,too-many-branches,too-many-locals,too-many-statements
        """
        WIG to HDF5 converter

        Loads the WIG file into the HDF5 index file that gets used by the REST
        API to determine if there are files that have data in a given region.
        Overlapping regions are condensed into a single feature block rather
        than maintaining all of the detail of the original WIG file.

        Parameters
        ----------
        file_id : str
            The file_id as stored by the DMP so that it can be used for file
            retrieval later
        assembly : str
            Assembly of the genome that is getting indexed so that the
            chromosomes match
        file_wig : str
            Location of the wig file
        file_hdf5 : str
            Location of the HDF5 index file

        Example
        -------
        .. code-block:: python
           :linenos:

           if not self.wig2hdf5(file_id, assembly, wig_file, hdf5_file):
               output_metadata.set_exception(
                   Exception(
                       "wig2hdf5: Could not process files {}, {}.".format(*input_files)))

        """
        max_files = 1024
        max_chromosomes = 1024
        max_chromosome_size = 2000000000

        hdf5_in = h5py.File(file_hdf5, "a")

        if str(assembly) in hdf5_in:
            grp = hdf5_in[str(assembly)]
            meta = hdf5_in['meta']

            dset = grp['data']
            fset = grp['files']
            cset = grp['chromosomes']
            file_idx = [i for i in fset if i != '']
            if file_id not in file_idx:
                file_idx.append(file_id)
                # pylint is unable to recognise the resize and shape methods
                dset.resize((dset.shape[0], dset.shape[1] + 1, max_chromosome_size))  # pylint: disable=no-member
            chrom_idx = [c for c in cset if c != '']

        else:
            # Create the initial dataset with minimum values
            grp = hdf5_in.create_group(str(assembly))
            meta = hdf5_in.create_group('meta')

            dtf = h5py.special_dtype(vlen=str)
            dtc = h5py.special_dtype(vlen=str)
            fset = grp.create_dataset('files', (max_files,), dtype=dtf)
            cset = grp.create_dataset('chromosomes', (max_chromosomes,), dtype=dtc)

            file_idx = [file_id]
            chrom_idx = []

            dset = grp.create_dataset(
                'data1', (0, 1, max_chromosome_size),
                maxshape=(max_chromosomes, max_files, max_chromosome_size),
                dtype='bool', chunks=True, compression="gzip")

        # Save the list of files
        fset[0:len(file_idx)] = file_idx

        file_chrom_count = 0

        dnp = np.zeros([max_chromosome_size], dtype='bool')

        loaded = False

        start = 0
        span = 1
        step = 1
        step_type = 'variable'
        previous_chrom = ''
        previous_start = 0
        previous_end = 0

        with open(file_wig, 'r') as f_in:
            for line in f_in:
                line = line.strip()
                if line[0:9] == 'fixedStep' or line[0:12] == 'variableStep':
                    start = 0
                    span = 1
                    step = 1
                    previous_chrom = ''
                    step_type = 'variable'

                    sline = line.split(' ')

                    if sline[0][0:9] == 'fixedStep':
                        step_type = 'fixed'

                    if previous_chrom != '':
                        dnp[previous_start:previous_end + 1] = 1

                    chrom = ''
                    for key_value in sline[1:]:
                        key_value.strip()
                        if not key_value:
                            continue
                        k, i = key_value.split('=')
                        if k == 'start':
                            start = int(i)
                        elif k == 'span':
                            span = int(i)
                        elif k == 'step':
                            step = int(i)
                        elif k == 'chrom':
                            chrom = i

                    file_chrom_count += 1
                    if previous_chrom != '' and previous_chrom != chrom:
                        if previous_chrom not in chrom_idx:
                            chrom_idx.append(previous_chrom)
                            cset[0:len(chrom_idx)] = chrom_idx
                            dset.resize((dset.shape[0] + 1, dset.shape[1], max_chromosome_size))

                        dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] += dnp

                        dnp = np.zeros([max_chromosome_size], dtype='bool')
                        loaded = True

                    previous_chrom = chrom

                else:
                    loaded = False
                    if step_type == 'fixed':
                        if float(line) == 0.0:
                            if previous_start != previous_end:
                                dnp[previous_start:previous_end + 1] = 1
                                previous_start = 0
                                previous_end = 0
                        else:
                            if previous_start == 0:
                                previous_start = start
                                previous_end = start + span - 1
                            else:
                                if previous_end == start - 1:
                                    previous_end += span
                                else:
                                    dnp[previous_start:previous_end + 1] = 1
                                    previous_start = start
                                    previous_end = start + span - 1

                        start += step

                    elif step_type == 'variable':
                        sline = line.split("\t")
                        if float(sline[1]) == 0.0:
                            if previous_start != previous_end:
                                dnp[previous_start:previous_end + 1] = 1
                                previous_start = 0
                                previous_end = 0
                        else:
                            if previous_start == 0:
                                previous_start = int(sline[0])
                                previous_end = int(sline[0]) + span - 1
                            else:
                                if previous_end == int(sline[0]) - 1:
                                    previous_end += span
                                else:
                                    dnp[previous_start:previous_end + 1] = 1
                                    previous_start = int(sline[0])
                                    previous_end = int(sline[0]) + span - 1

            if loaded is False:
                if previous_chrom not in chrom_idx:
                    chrom_idx.append(chrom)
                    cset[0:len(chrom_idx)] = chrom_idx
                    dset.resize((dset.shape[0] + 1, dset.shape[1], max_chromosome_size))

                dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp

        hdf5_in.close()

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Function to run the WIG file sorter and indexer so that the files can
        get searched as part of the REST API

        Parameters
        ----------
        input_files : dict
            wig_file : str
                Location of the wig file
            chrom_size : str
                Location of chrom.size file
            hdf5_file : str
                Location of the HDF5 index file
        meta_data : dict

        Returns
        -------
        list
            bw_file : str
                Location of the BigWig file
            hdf5_file : str
                Location of the HDF5 index file

        """
        logger.info(
            'PIPELINE FILES:', input_files["wig"], input_files["chrom_file"],
            output_files["bw_file"])
        results_1 = self.wig2bigwig(
            input_files["wig"], input_files["chrom_file"], output_files["bw_file"])
        results_1 = compss_wait_on(results_1)

        results_2 = self.wig2hdf5(
            input_files["wig"], input_metadata["wig"].meta_data["assembly"],
            input_files["wig"], input_files["hdf5_file"])
        results_2 = compss_wait_on(results_2)

        output_generated_files = {
            "bw_file": output_files["bw_file"],
            "hdf5_file": input_files["hdf5_file"]
        }
        output_metadata = {
            "bw_file": input_files["wig"],
            "hdf5_file": input_metadata["hdf5_file"]
        }

        return (output_generated_files, output_metadata)

# ------------------------------------------------------------------------------
