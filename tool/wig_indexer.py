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

class wigIndexerTool(Tool):
    """
    Tool for running indexers over a WIG file for use in the RESTful API
    """
    
    def __init__(self):
        """
        Init function
        """
        print "WIG File Indexer"
    
    
    @task(file_wig=FILE_IN, file_chrom=FILE_IN, file_bw=FILE_OUT)
    def wig2bigwig(self, file_wig, file_chrom, file_bw):
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
        command_line = 'wigToBigWig ' + file_wig + ' ' + file_chrom + ' ' + file_bw
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        return True
    
    
    @task(file_id=IN, assembly=IN, file_wig=FILE_IN, file_hdf5=FILE_INOUT)
    def wig2hdf5(self, file_id, assembly, file_wig, file_hdf5):
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
        MAX_FILES = 1024
        MAX_CHROMOSOMES = 1024
        MAX_CHROMOSOME_SIZE = 2000000000
        
        f = h5py.File(resource_path, "a")
        
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

        loaded = False
        
        fi = open(file_wig, 'r')
        start = 0
        span  = 1
        step  = 1
        step_type = 'variable'
        previous_chrom = ''
        previous_start = 0
        previous_end   = 0

        for line in fi:
            line = line.strip()
            #print(start)
            if line[0:9] == 'fixedStep' or line[0:12] =='variableStep':
                start = 0
                span  = 1
                step  = 1
                previous_chrom = ''
                step_type = 'variable'
                
                l = line.split(' ')
                
                if l[0][0:9] == 'fixedStep':
                    step_type = 'fixed'
                
                if previous_chrom != '':
                    dnp[previous_start:previous_end+1] = 1
                
                chrom = ''
                for kv in l[1:]:
                    kv.strip()
                    if len(kv) == 0:
                        continue
                    k, v = kv.split('=')
                    if k == 'start':
                        start = int(v)
                    elif k == 'span':
                        span = int(v)
                    elif k == 'step':
                        step = int(v)
                    elif k == 'chrom':
                        chrom = v
                
                file_chrom_count += 1
                if previous_chrom != '' and previous_chrom != chrom:
                    if previous_chrom not in chrom_idx:
                        chrom_idx.append(previous_chrom)
                        cset[0:len(chrom_idx)] = chrom_idx
                        dset.resize((dset.shape[0]+1, dset.shape[1], MAX_CHROMOSOME_SIZE))
                    
                    dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] += dnp
                    
                    dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='bool')
                    loaded = True
                
                previous_chrom = chrom
                
            else:
                loaded = False
                #print(chrom, str(start), str(step), str(span))
                if step_type == 'fixed':
                    if float(line) == 0.0:
                        if previous_start != previous_end:
                            dnp[previous_start:previous_end+1] = 1
                            previous_start = 0
                            previous_end   = 0
                    else:
                        if previous_start == 0:
                            previous_start = start
                            previous_end   = start + span - 1
                        else:
                            if previous_end == start-1:
                                previous_end += span
                            else:
                                dnp[previous_start:previous_end+1] = 1
                                previous_start = start
                                previous_end   = start + span - 1
                    
                    start += step
                    
                elif step_type == 'variable':
                    l = line.split("\t")
                    if float(l[1]) == 0.0:
                        if previous_start != previous_end:
                            dnp[previous_start:previous_end+1] = 1
                            previous_start = 0
                            previous_end   = 0
                    else:
                        if previous_start == 0:
                            previous_start = int(l[0])
                            previous_end   = int(l[0]) + span - 1
                        else:
                            if previous_end == int(l[0])-1:
                                previous_end += span
                            else:
                                dnp[previous_start:previous_end+1] = 1
                                previous_start = int(l[0])
                                previous_end   = int(l[0]) + span - 1
        
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
        Function to run the WIG file sorter and indexer so that the files can
        get searched as part of the REST API
        
        Parameters
        ----------
        input_files : list
            wig_file : str
                Location of the wig file
            chrom_size : str
                Location of chrom.size file
            hdf5_file : str
                Location of the HDF5 index file
        meta_data : list
            file_id : str
                file_id used to identify the original wig file
            assembly : str
                Genome assembly accession
        
        Returns
        -------
        list
            bw_file : str
                Location of the BigWig file
            hdf5_file : str
                Location of the HDF5 index file
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           import tool
           
           # WIG Indexer
           w = tool.bedIndexerTool(self.configuration)
           wi, wm = w.run((wig_file_id, chrom_file_id, hdf5_file_id), {'file_id' : file_id, 'assembly' : assembly})
        """
        wig_file   = input_files[0]
        chrom_file = input_files[1]
        hdf5_file  = input_files[2]
        
        bw_name = wig_file.split("/")
        bw_name[-1].replace('.wig', '.bw')
        bw_file = '/'.join(bw_name)
        
        assembly = meta_data['assembly']
        
        # handle error
        if not self.wig2bigwig(wig_file, chrom_file, bw_file):
            output_metadata.set_exception(
                Exception(
                    "wig2bigwig: Could not process files {}, {}.".format(*input_files)))
        
        if not self.wig2hdf5(file_id, assembly, wig_file, hdf5_file):
            output_metadata.set_exception(
                Exception(
                    "wig2hdf5: Could not process files {}, {}.".format(*input_files)))
        
        return ([bw_file, hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
