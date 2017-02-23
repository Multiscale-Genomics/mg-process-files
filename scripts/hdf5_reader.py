#!/usr/bin/python

"""
Copyright 2017 EMBL-European Bioinformatics Institute

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
    
    
    def __del__(self):
        """
        Tidy function to close file handles
        """
        self.f.close()
    
    
    def close(self):
        """
        Tidy function to close file handles
        """
        self.f.close()
    
    
    def get_assemblies(self):
        """
        List all assemblies for which there are files that have been indexed
        """
        return [asm for asm in self.f[self.user_id] if asm != 'meta']
    
    
    def get_chromosomes(self):
        """
        List all chromosomes that are covered by the index
        """
        grp = self.f[str(self.user_id)]
        cset = grp['chromosomes']
        return [c for c in cset if c != '']
    
    
    def get_files(self):
        """
        List all files for an assembly. If files are missing they can either get
        loaded or the search can be performed directly on the bigBed files
        """
        grp = self.f[str(self.user_id)]
        fset = grp['files']
        return [f for f in fset if f != '']
    
    
    def get_regions(self, assembly, chromosome_id, start, end):
        """
        List files that have data in a given region.
        """
        grp = self.f[str(self.user_id)][assembly]
        dset = grp['data']
        fset = grp['files']
        cset = grp['chromosomes']
        
        fid = list(np.nonzero(fset))
        file_idx = [ cset[i] for i in fid[0] ]
        
        cid = list(np.nonzero(cset))
        chrom_idx = [ cset[i] for i in cid[0] ]
        
        c = str(chromosome_id)
        s = int(start)
        e = int(end)
        
        #cfp
        dnp = dset[chrom_idx.index(c),:,s:e]
        f_idx = []
        for i in range(len(dnp)):
            if np.any(dnp[i]) == True:
                f_idx.append(file_idx[i])
        return f_idx
    
