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

import pytadbit
from pytadbit.parsers.hic_parser import read_matrix
import numpy as np
import h5py

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT, FILE_INOUT
    from pycompss.api.task import task
except ImportError :
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")
    
    from dummy_pycompss import *

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class hicAdjacencyIndexerTool(Tool):
    """
    Tool for running indexers over a BED file for use in the RESTful API
    """
    
    def __init__(self):
        """
        Init function
        """
        print "Adjacency Indexer"
    

    @task(adjlist_file = FILE_IN, hdf5_file = FILE_INOUT)
    def hic2hdf5(self, adjlist_file, hdf5_file, resolution):
        """
        Save adjacency matrix as HDF5 file

        Save the hic_data object to HDF5 file. This is saved as an NxN array
        with the values for all positions being set.
        
        This needs to include attributes for the chromosomes for each resolution
        - See the mg-rest-adjacency hdf5_reader for further details about the
          requirement. This prevents the need for secondary storage details
          outside of the HDF5 file.

        Parameters
        ----------
        adjlist_file : str
            Location of the adjacency list file
        hdf5_file : str
            Location of the HDF5 index file
        resolution : int
            Resolution of the adjacency list calls
        """
        hic_data = read_matrix(adjlist_file, resolution=resolution)

        dSize = len(hic_data)
        d = np.zeros([dSize, dSize], dtype='int32')
        d += hic_data.get_matrix()
        
        f = h5py.File(hdf5_file, "a")
        dset = f.create_dataset(str(self.resolution), (dSize, dSize), dtype='int32', chunks=True, compression="gzip")
        dset[0:dSize,0:dSize] += d
        f.close()

        return True


    def run(self, input_files, metadata):
        """
        Function to store the TADbit generated adjacency list data within an
        HDF5 index for use from the REST interface
        
        Parameters
        ----------
        input_files : list
            adj_list_file : str
                Location of the adjacency list file
            hdf5_file : str
                Location of the HDF5 index file
        meta_data : list
            file_id : str
                file_id used to identify the original bed file
        
        Returns
        -------
        list
            hdf5_file : str
                Location of the HDF5 index file
        
        Example
        -------
        .. code-block:: python
           :linenos:
           
           import tool
           
           # Adjacency File Indexer
           al = tool.hicAdjacencyIndexerTool(self.configuration)
           ali, alm = al.run((adjlist_file, hdf5_file), {'assembly' : assembly, 'resolution' : resolution})
        """
        adjlist_file = input_files[0]
        hdf5_file    = input_files[1]
        
        assembly = metadata['assembly']
        resolution = metadata['resolution']
        
        # handle error
        if not self.hic2hdf5(adjlist_file, hdf5_file, resolution):
            output_metadata.set_exception(
                Exception(
                    "hic2hdf5: Could not process files {}, {}.".format(*input_files)))
        
        return ([hdf5_file], [output_metadata])

# ------------------------------------------------------------------------------
