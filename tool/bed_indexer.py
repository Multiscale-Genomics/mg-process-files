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

class bedIndexerTool(Tool):
    """
    Tool for running indexers over a BED file for use in the RESTful API
    """
    
    def __init__(self):
        """
        Init function
        """
        print "BED File Indexer"

        self.feature_break_length = 50
    
    
    @task(file_bed=FILE_IN, file_sorted_bb=FILE_OUT)
    def bedsort(self, file_bed, file_sorted_bed):
        """
        BED file sorter
        
        This is a wrapper for the standard Linux ``sort`` method the sorting by
        the chromosome and start columns in the BED file.
        
        Parameters
        ----------
        file_bed : str
            Location of the BED file
        file_sorted_bed : str
            Location of the sorted BED file
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           if not self.bedsorted(bed_file, bed_sorted_file):
               output_metadata.set_exception(
                   Exception(
                       "bedsorted: Could not process files {}, {}.".format(*input_files)))
        
        """
        with open(file_sorted_bed,"wb") as out:
            command_line = 'sort -k1,1 -k2,2n ' + file_bed
            args = shlex.split(command_line)
            p = subprocess.Popen(args,stdout=out)
            p.wait()
        return True
    
    
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

        # Features with length <100bp
        small_feature_count = 0
        small_feature_length = 0

        # Features with length >= 100bp
        long_feature_count = 0
        long_feature_length = 0

        fi = open(file_bed, 'r')
        for line in fi:
            line = line.strip()
            l = line.split("\t")
            
            c = str(l[0])
            s = int(l[1])
            e = int(l[2])
            l = e-s

            if l < self.feature_break_length:
                small_feature_count += 1
                small_feature_length += l
            else:
                long_feature_count += 1
                long_feature_length += l

        total_feature_length = small_feature_length + long_feature_length
        total_feature_count = small_feature_count + long_feature_count
        return total_feature_length / total_feature_count

    @task(file_sorted_bed=FILE_IN, file_chrom=FILE_IN, file_bb=FILE_OUT, bed_type=IN)
    def bed2bigbed(self, file_sorted_bed, file_chrom, file_bb, bed_type = None):
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
        if bed_type != None:
            command_line += ' -type=' + str(bed_type)
        command_line +=  ' ' + file_sorted_bed + ' ' + file_chrom + ' ' + file_bb
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        return True
    
    
    @task(file_id=IN, assembly=IN, feature_length=IN, file_sorted_bed=FILE_IN, file_hdf5=FILE_INOUT)
    def bed2hdf5(self, file_id, assembly, feature_length, file_sorted_bed, file_hdf5):
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
        MAX_FILES = 1024
        MAX_CHROMOSOMES = 1024
        MAX_CHROMOSOME_SIZE = 2000000000
        
        f = h5py.File(file_hdf5, "a")
        
        if str(assembly) in f:
            grp  = f[str(assembly)]
            meta = f['meta']
            
            dset1  = grp['data1']
            dset1k = grp['data1k']
            fset = grp['files']
            cset = grp['chromosomes']
            file_idx_1  = [fs for fs in fset[0] if fs != '']
            file_idx_1k = [fs for fs in fset[1] if fs != '']
            if file_id not in file_idx_1 and file_id not in file_idx_1k:
                if feature_length == 1000:
                    file_idx_1k.append(file_id)
                else:
                    file_idx_1.append(file_id)
                dset1.resize((dset1.shape[0], dset1.shape[1]+1, MAX_CHROMOSOME_SIZE))
                dset1k.resize((dset1k.shape[0], dset1k.shape[1]+1, MAX_CHROMOSOME_SIZE/1000))
            chrom_idx = [c for c in cset if c != '']
                
        else:
            # Create the initial dataset with minimum values
            grp    = f.create_group(str(assembly))
            meta   = f.create_group('meta')
            
            dtf = h5py.special_dtype(vlen=str)
            dtc = h5py.special_dtype(vlen=str)
            fset = grp.create_dataset('files', (2, MAX_FILES), dtype=dtf)
            cset = grp.create_dataset('chromosomes', (MAX_CHROMOSOMES,), dtype=dtc)
            
            file_idx_1  = []
            file_idx_1k = []
            chrom_idx = []
            
            dset1 = grp.create_dataset('data1', (0, 1, MAX_CHROMOSOME_SIZE),
                maxshape=(MAX_CHROMOSOMES, MAX_FILES, MAX_CHROMOSOME_SIZE),
                dtype='str', chunks=True, compression="gzip")
            dset1k = grp.create_dataset('data1k', (0, 1, MAX_CHROMOSOME_SIZE/1000),
                maxshape=(MAX_CHROMOSOMES, MAX_FILES, MAX_CHROMOSOME_SIZE/1000),
                dtype='str', chunks=True, compression="gzip")

            if feature_length == 1000:
                file_idx_1k.append(file_id)
            else:
                file_idx_1.append(file_id)
        
        # Save the list of files
        fset[0, 0:len(file_idx_1)] = file_idx_1
        fset[0, 0:len(file_idx_1k)] = file_idx_1k
        
        file_chrom_count = 0

        if feature_length == 1000:
            # dnp = np.full([MAX_CHROMOSOME_SIZE/1000], None, dtype='str')
            dnp = np.zeros([MAX_CHROMOSOME_SIZE/1000], dtype='int8')
        else:
            # dnp = np.full([MAX_CHROMOSOME_SIZE], None, dtype='str')
            dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='int8')

        previous_chrom = ''
        previous_start = 0
        previous_end   = 0
        
        loaded = False
        
        fi = open(file_sorted_bed, 'r')
        for line in fi:
            line = line.strip()
            l = line.split("\t")
            
            c = str(l[0])
            s = int(l[1])
            e = int(l[2])
            
            loaded = False
            
            if c != previous_chrom and previous_chrom != '':
                #if c != 'chr2':
                #    loaded = True
                #    continue
                file_chrom_count += 1
                if previous_chrom not in chrom_idx:
                    chrom_idx.append(previous_chrom)
                    cset[0:len(chrom_idx)] = chrom_idx
                    dset1.resize((dset1.shape[0]+1, dset1.shape[1], MAX_CHROMOSOME_SIZE))
                    dset1k.resize((dset1k.shape[0]+1, dset1k.shape[1], MAX_CHROMOSOME_SIZE/1000))
                
                loaded = True
                
                if feature_length == 1000:
                    dset1k[chrom_idx.index(previous_chrom), file_idx_1k.index(file_id), :] = dnp
                    # dnp = np.full([MAX_CHROMOSOME_SIZE/1000], None, dtype='str')
                    dnp = np.zeros([MAX_CHROMOSOME_SIZE/1000], dtype='int8')
                else:
                    dset1[chrom_idx.index(previous_chrom), file_idx_1.index(file_id), :] = dnp
                    # dnp = np.full([MAX_CHROMOSOME_SIZE], None, dtype='str')
                    dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='int8')
            
            previous_chrom = c
            if feature_length == 1000:
                dnp[(s/1000):(e/1000)+1] = '1'
            else:
                dnp[s:e+1] = '1'

        if loaded == False:
            if previous_chrom not in chrom_idx:
                chrom_idx.append(c)
                cset[0:len(chrom_idx)] = chrom_idx
                dset1.resize((dset1.shape[0]+1, dset1.shape[1], MAX_CHROMOSOME_SIZE))
                dset1k.resize((dset1k.shape[0]+1, dset1k.shape[1], MAX_CHROMOSOME_SIZE/1000))
            
            if feature_length == 1000:
                dset1k[chrom_idx.index(previous_chrom), file_idx_1k.index(file_id), :] = dnp
            else:
                dset1[chrom_idx.index(previous_chrom), file_idx_1.index(file_id), :] = dnp
        
        f.close()
        fi.close()
        
        return True
    
    
    def run(self, input_files, metadata):
        """
        Function to run the BED file sorter and indexer so that the files can
        get searched as part of the REST API
        
        Parameters
        ----------
        input_files : list
            bed_file : str
                Location of the bed file
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
           bi, bm = bd.run((bed_file_id, chrom_file_id, hdf5_file_id), {'file_id' : file_id, 'assembly' : assembly})
        """
        bed_file   = input_files[0]
        chrom_file = input_files[1]
        hdf5_file  = input_files[2]
        
        bed_sorted_name = bed_file.split("/")
        bed_sorted_name[-1] = bed_sorted_name[-1].replace('.bed', '.sorted.bed')
        bed_sorted_file = '/'.join(bed_sorted_name)
        
        bb_name = bed_file.split("/")
        bb_name[-1] = bb_name[-1].replace('.bed', '.bb')
        bb_file = '/'.join(bb_name)
        
        assembly = metadata['assembly']
        bed_type = metadata['bed_type']

        output_metadata = {}

        # handle error
        if not self.bedsort(bed_file, bed_sorted_file):
            output_metadata.set_exception(
                Exception(
                    "bedsorted: Could not process files {}, {}.".format(*input_files)))
        
        if not self.bed2bigbed(bed_file, chrom_file, bb_file, bed_type):
            output_metadata.set_exception(
                Exception(
                    "bed2bigbed: Could not process files {}, {}.".format(*input_files)))
        
        feature_length = self.bed_feature_length(bed_sorted_file)
        storage_level = 1000
        if feature_length < 10:
            storage_level = 1

        if not self.bed2hdf5(metadata['file_id'], assembly, storage_level, bed_sorted_file, hdf5_file):
            output_metadata.set_exception(
                Exception(
                    "bed2hdf5: Could not process files {}, {}.".format(*input_files)))
        
        return ([bb_file, hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
