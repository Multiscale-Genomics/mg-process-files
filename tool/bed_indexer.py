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
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
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
    
    
    @task(file_sorted_bed=FILE_IN, file_chrom=FILE_IN, file_bb=FILE_OUT)
    def bed2bigbed(self, file_bed, file_chrom, file_bb):
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
        command_line = 'bedToBigBed ' + file_sorted_bed + ' ' + file_chrom + ' ' + file_bb
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        return True
    
    
    @task(file_id=IN, assembly=IN, file_sorted_bed=FILE_IN, file_hdf5=FILE_INOUT)
    def bed2hdf5(self, file_id, assembly, file_sorted_bed, file_hdf5):
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
        
        fi = open(file_sorted_bed, 'r')
        for line in fi:
            line = line.strip()
            l = line.split("\t")
            
            c = str(l[0])
            s = int(l[1])
            e = int(l[2])
            
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
        meta_data : list
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
        bed_sorted_name[-1].replace('.bed', '.sorted.bed')
        bed_sorted_file = '/'.join(bed_sorted_name)
        
        bb_name = bed_file.split("/")
        bb_name[-1].replace('.bed', '.bb')
        bb_file = '/'.join(bb_name)
        
        assembly = meta_data['assembly']
        
        # handle error
        if not self.bedsorted(bed_file, bed_sorted_file):
            output_metadata.set_exception(
                Exception(
                    "bedsorted: Could not process files {}, {}.".format(*input_files)))
        
        if not self.bed2bigbed(bed_file, chrom_file, bb_file):
            output_metadata.set_exception(
                Exception(
                    "bed2bigbed: Could not process files {}, {}.".format(*input_files)))
        
        if not self.bed2hdf5(file_id, assembly, bed_file, hdf5_file):
            output_metadata.set_exception(
                Exception(
                    "bed2hdf5: Could not process files {}, {}.".format(*input_files)))
        
        return ([bb_file, hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
