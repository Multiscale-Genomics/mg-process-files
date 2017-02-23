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

import argparse, urllib2, gzip, shutil, shlex, subprocess, os.path, json

from basic_modules import Tool, Workflow, Metadata

from dmp import dmp

import tool

import os

try
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.constraint import constraint
except ImportError :
    print "[Warning] Cannot import \"pycompss\" API packages."
    print "          Using mock decorators."
    
    from dummy_pycompss import *

try :
    import pysam
except ImportError :
    print "[Error] Cannot import \"pysam\" package. Have you installed it?"
    exit(-1)

# ------------------------------------------------------------------------------

class processs_bed(Workflow):
    """
    Workflow to index bed formatted files
    """
    
    def run(self, file_ids, metadata):
        """
        Main run function
        
        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list
        
        Returns
        -------
        outputfiles : list
            List of locations for the output bam, bed and tsv files
        """
        
        bed_fa = file_ids[0]
        chrom_file = file_ids[1]
        hdf5_file = file_ids[2]
        assembly = metadata["assembly"]
        
        # Bed Indexer
        b = tool.bedIndexerTool(self.configuration)
        bb, h5_idx = bd.run((bed_file, chrom_file, hdf5_file), {'assembly' : assembly})
        
        return (bb, h5_idx)

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os
    
    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Index the bed file")
    parser.add_argument("--assembly", help="Assembly")
    parser.add_argument("--chrom", help="Matching chrom.size file")
    parser.add_argument("--bed_file", help="Bed file to get indexed")
    parser.add_argument("--h5_file", help="Bed file to get indexed")
    
    # Get the matching parameters from the command line
    args = parser.parse_args()
    
    assembly = args.assembly
    chrom_size_file = args.chrom
    bed_file = args.bed_file
    h5_file = args.h5_file
    
    pb = process_bed()
    
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
    
    cs_file = da.set_file("test", chrom_size_file, "tsv", "ChIP-seq", "", None)
    b_file = da.set_file("test", bed_file, "bed", "Assembly", "", None)
    h5_file = da.set_file("test", h5_file, "hdf5", "ChIP-seq", "", None)
    
    print da.get_files_by_user("test")
    
    # 3. Instantiate and launch the App
    from basic_modules import WorkflowApp
    app = WorkflowApp()
    results = app.launch(process_bed, [b_file, cs_file, h5_file], {"assembly" : assembly})
    
    print da.get_files_by_user("test")
