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

import argparse, gzip, shutil, shlex, subprocess, os.path, json
from functools import wraps


from basic_modules import Tool, Workflow, Metadata
from dmp import dmp

import tool
import os

try:
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
except ImportError :
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")
    
    from dummy_pycompss import *


# ------------------------------------------------------------------------------

class process_adjacency_list(Workflow):
    """
    Workflow to index BED formatted files within the Multiscale Genomics (MuG)
    Virtural Research Environment (VRE)
    """
    
    def __init__(self):
        """
        Initialise the class
        """
        
    
    def run(self, file_ids, metadata):
        """
        Main run function to index the BED files ready for use in the RESTful
        API. BED files are index in 2 different ways to allow for optimal data
        retreival. The first is as a bigbed file, this allows the data to get
        easily extracted as BED documents and served to the user. The second is
        as an HDF5 file that is used to identify which bed files have
        information at a given location. This is to help the REST clients make
        only the required calls to the relevant BED files rather than needing to
        pole all potential BED files.
        
        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
        
        Returns
        -------
        outputfiles : list
            List of locations for the output BED and HDF5 files
        """
        
        bed_file = file_ids[0]
        chrom_file = file_ids[1]
        hdf5_file = file_ids[2]
        assembly = metadata["assembly"]
        
        # Bed Indexer
        b = tool.bedIndexerTool(self.configuration)
        bb, h5_idx = b.run((bed_file, chrom_file, hdf5_file), {'assembly' : assembly})
        
        return (bb, h5_idx)

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Index the TADbit adjacency file")
    parser.add_argument("--assembly", help="Assembly")
    parser.add_argument("--adjlist_file", help="Adjacency list file to get indexed")
    parser.add_argument("--h5_file", help="Location of HDF5 index file")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    assembly = args.assembly
    adjlist_file = args.adjlist_file
    h5_file = args.h5_file
    
    #
    # MuG Tool Steps
    # --------------
    # 
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    Are ones that are present in the ENA then I would need to download them
    
    
    #2. Register the data with the DMP
    da = dmp()
    
    print da.get_files_by_user("test")
    
    a_file = da.set_file("test", adjlist_file, "tsv", "Hi-C", "", {'assembly' : assembly})
    h5_file = da.set_file("test", h5_file, "hdf5", "index", "", None)
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    from basic_modules import WorkflowApp
    app = WorkflowApp()
    results = app.launch(process_adjacency_list, [a_file, h5_file], {"assembly" : assembly})
    
    print da.get_files_by_user("test")
