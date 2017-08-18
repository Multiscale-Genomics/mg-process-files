"""
.. Copyright 2017 EMBL-European Bioinformatics Institute

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

import os
import subprocess
import shlex

import numpy as np
import h5py

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
except ImportError :
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class gff3IndexerTool(Tool):
    """
    Tool for running indexers over a WIG file for use in the RESTful API
    """

    def __init__(self):
        """
        Init function
        """
        print "GFF3 File Indexer"


    def gff3Sorter(self, file_gff3):
        """
        Sorts the GFF3 file in preparation for compression and indexing

        Parameters
        ----------
        file_gff3 : str
            Location of the source GFF3 file
        """
        grep_comment_command = 'grep ^"#" ' + file_gff3
        grep_gff_command = 'grep -v ^"#" ' + file_gff3
        sort_command = 'sort -k1,1 -k4,4n'

        grep_comment_args = shlex.split(grep_comment_command)
        grep_args = shlex.split(grep_gff_command)
        sort_args = shlex.split(sort_command)

        file_sorted_gff3 = file_gff3 + '.sorted.gff3'

        with open(file_sorted_gff3,"wb") as out:
            grep = subprocess.Popen(grep_comment_args, stdout=out)
            grep.wait()

        with open(file_sorted_gff3,"a") as out:
            grep = subprocess.Popen(grep_args, stdout=subprocess.PIPE)
            sorted = subprocess.Popen(sort_args, stdin=grep.stdout, stdout=out)
            sorted.wait()

        return file_sorted_gff3


    @task(file_sorted_gff3=FILE_IN, file_sorted_gz_gff3=FILE_OUT, file_gff3_tbi=FILE_OUT)
    def gff32tabix(self, file_sorted_gff3, file_sorted_gz_gff3, file_gff3_tbi):
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
        pysam.tabix_compress(file_sorted_gff3, file_sorted_gz_gff3)
        pysam.tabix_index(file_sorted_gz_gff3, preset='gff')
        return True

    @task(file_id=IN, assembly=IN, file_sorted_gff3=FILE_IN, file_hdf5=FILE_INOUT)
    def gff32hdf5(self, file_id, assembly, file_sorted_gff3, file_hdf5):
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
        MAX_FILES = 1024
        MAX_CHROMOSOMES = 1024
        MAX_CHROMOSOME_SIZE = 2000000000

        f = h5py.File(file_hdf5, "a")

        if str(assembly) in f:
            grp  = f[str(assembly)]
            meta = f['meta']

            dset = grp['data']
            fset = grp['files']
            cset = grp['chromosomes']
            file_idx  = [f for f in fset if f != '']
            if file_id not in file_idx:
                file_idx.append(file_id)
                dset.resize((dset.shape[0], dset.shape[1]+1, MAX_CHROMOSOME_SIZE))
            chrom_idx = [c for c in cset if c != '']

        else:
            # Create the initial dataset with minimum values
            grp    = f.create_group(str(assembly))
            meta   = f.create_group('meta')

            dtf = h5py.special_dtype(vlen=str)
            dtc = h5py.special_dtype(vlen=str)
            fset = grp.create_dataset('files', (MAX_FILES,), dtype=dtf)
            cset = grp.create_dataset('chromosomes', (MAX_CHROMOSOMES,), dtype=dtc)

            file_idx  = [file_id]
            chrom_idx = []

            dset = grp.create_dataset('data', (0, 1, MAX_CHROMOSOME_SIZE), maxshape=(MAX_CHROMOSOMES, MAX_FILES, MAX_CHROMOSOME_SIZE), dtype='bool', chunks=True, compression="gzip")

        # Save the list of files
        fset[0:len(file_idx)] = file_idx

        file_chrom_count = 0

        dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='bool')

        previous_chrom = ''
        previous_start = 0
        previous_end   = 0

        loaded = False

        fi = open(file_sorted_gff3, 'r')
        for line in fi:
            line = line.strip()
            l = line.split("\t")

            c = str(l[0])
            s = int(l[3])
            e = int(l[4])

            loaded = False

            if c != previous_chrom and previous_chrom != '':
                file_chrom_count += 1
                if previous_chrom not in chrom_idx:
                    chrom_idx.append(previous_chrom)
                    cset[0:len(chrom_idx)] = chrom_idx
                    dset.resize((dset.shape[0]+1, dset.shape[1], MAX_CHROMOSOME_SIZE))

                dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp
                loaded = True

                if file_chrom_count == 5:
                    break

                dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='bool')

            previous_chrom = c
            dnp[s:e+1] = 1

        if loaded == False:
            if previous_chrom not in chrom_idx:
                chrom_idx.append(c)
                cset[0:len(chrom_idx)] = chrom_idx
                dset.resize((dset.shape[0]+1, dset.shape[1], MAX_CHROMOSOME_SIZE))

            dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp

        f.close()
        fi.close()

        return True


    def run(self, input_files, output_files, metadata=None):
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

        Example
        -------
        .. code-block:: python
           :linenos:

           import tool

           # Bed Indexer
           g = tool.gff3IndexerTool(self.configuration)
           gff3_files, gff3_meta = g.run((gff3_file_id, hdf5_file_id), {'file_id' : file_id, 'assembly' : assembly})
        """
        gff3_file   = input_files[0]
        hdf5_file  = input_files[2]

        gz_file = gff3_file + '.gz'
        tbi_file = gz_file + '.tbi'

        assembly = meta_data['assembly']

        file_sorted_gff3 = self.gff3Sorter(gff3_file)

        # handle error
        if not self.gff32tabix(self, file_sorted_gff3, gz_file, tbi_file):
            output_metadata.set_exception(
                Exception(
                    "gff32tabix: Could not process files {}, {}.".format(*input_files)))
            gz_file  = None
            tbi_file = None

        if not self.gff32hdf5(file_id, assembly, file_sorted_gff3, hdf5_file):
            output_metadata.set_exception(
                Exception(
                    "gff32hdf5: Could not process files {}, {}.".format(*input_files)))

        return ([gz_file, tbi_file, hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
