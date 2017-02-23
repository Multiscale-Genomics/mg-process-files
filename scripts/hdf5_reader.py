#!/usr/bin/python

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

import os, json, h5py
import numpy as np

from dmp import dmp

class hdf5_reader:
    """
    Class related to handling the functions for interacting directly with the
    HDF5 files. All required information should be passed to this class.
    """
    
    test_file = '../region_idx.hdf5'
    user_id = 'test_user'
    
    def __init__(self, user_id = 'test'):
        """
        Initialise the module and 
        
        Parameters
        ----------
        user_id : str
            Identifier to uniquely locate the users files. Can be set to 
            "common" if the files can be shared between users
        file_id : str
            Location of the file in the file system
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           from hdf5_reader import hdf5_reader
           h5r = hdf5_reader('test')
        """
        
        self.test_file = 'region_idx.hdf5'
        self.user_id = 'test_user'
        
        # Open the hdf5 file
        if user_id == 'test':
            self.user_id = 'test_user'
            resource_path = os.path.join(os.path.dirname(__file__), self.test_file)
            self.f = h5py.File(resource_path, "r")
        else:
            self.user_id = user_id
            cnf_loc=os.path.dirname(os.path.abspath(__file__)) + '/mongodb.cnf'
            da = dmp(cnf_loc)
            file_obj = da.get_file_by_id(user_id, file_id)
            self.f = h5py.File(file_obj['file_path'], 'r')
    
    
    def close(self):
        """
        Tidy function to close file handles
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           from hdf5_reader import hdf5_reader
           h5r = hdf5_reader('test')
           h5r.close()
        """
        self.f.close()
    
    
    def get_assemblies(self):
        """
        List all assemblies for which there are files that have been indexed
        
        Returns
        -------
        assembly : list
            List of assemblies in the index
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           from hdf5_reader import hdf5_reader
           h5r = hdf5_reader('test')
           h5r.assemblies()
        """
        return [asm for asm in self.f[self.user_id] if asm != 'meta']
    
    
    def get_chromosomes(self, assembly):
        """
        List all chromosomes that are covered by the index
        
        Parameters
        ----------
        assembly : str
            Genome assembly ID
        
        Returns
        -------
        chromoosomes : list
            List of the chromosomes for a given assembly in the index
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           from hdf5_reader import hdf5_reader
           h5r = hdf5_reader('test')
           asm = h5r.assemblies()
           chr_list = h5r.get_chromosomes(asm[0])
        """
        grp = self.f[str(self.user_id)][assembly]
        cid = list(np.nonzero(grp['chromosomes']))
        return [ cset[i] for i in cid[0] ]
    
    
    def get_files(self, assembly):
        """
        List all files for an assembly. If files are missing they can either get
        loaded or the search can be performed directly on the bigBed files
        
        Parameters
        ----------
        assembly : str
            Genome assembly ID
        
        Returns
        -------
        file_ids : list
            List of file ids for a given assembly in the index
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           from hdf5_reader import hdf5_reader
           h5r = hdf5_reader('test')
           asm = h5r.assemblies()
           file_list = h5r.get_files(asm[0])
        """
        grp = self.f[str(self.user_id)][assembly]
        fid = list(np.nonzero(grp['files']))
        return [ fset[i] for i in fid[0] ]
    
    
    def get_regions(self, assembly, chromosome_id, start, end):
        """
        List files that have data in a given region.
        
        Parameters
        ----------
        assembly : str
            Genome assembly ID
        chromosome_id : str
            Chromosome names as listed by the get_files function
        start : int
            Start position for the region of interest
        end : int
            End position for the region of interest
        
        Returns
        -------
        file_ids : list
            List of the file_ids that have sequence features within the region
            of interest
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           from hdf5_reader import hdf5_reader
           h5r = hdf5_reader('test')
           asm = h5r.assemblies()
           file_list = h5r.get_chromosomes(asm[0], 1, 1000000, 1100000)
        """
        file_idx = self.get_files(assembly)
        chrom_idx = self.get_chromosomes(assembly)
        
        grp = self.f[str(self.user_id)][assembly]
        dset = grp['data']
        
        c = str(chromosome_id)
        s = int(start)
        e = int(end)
        
        #chr X file X position
        dnp = dset[chrom_idx.index(c),:,s:e]
        f_idx = []
        for i in range(len(dnp)):
            if np.any(dnp[i]) == True:
                f_idx.append(file_idx[i])
        return f_idx
    
