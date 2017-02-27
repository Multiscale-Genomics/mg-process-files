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

class json3dIndexerTool(Tool):
    """
    Tool for running indexers over 3D JSON files for use in the RESTful API
    """
    
    def __init__(self):
        """
        Init function
        """
        print "WIG File Indexer"
    
    
    
    def unzipJSON(self, file_targz):
        """
        
        """
        command_line = 'tar -xzf ' + file_tar
        args = shlex.split(command_line)
        p = subprocess.Popen(args)
        p.wait()
        return True
    
    
    def json2hdf5(self, source_file_id, targz_file, hdf5_file):
        """
        
        """
        
    
    
    def run(self, input_files, metadata):
        """
        
        """
        
        targz_file   = input_files[0]
        
        hdf5_name = targz_file.split("/")
        hdf5_name[-1].replace('.tar.gz', '.hdf5')
        hdf5_file = '/'.join(hdf5_name)
        
        source_file_id = meta_data['file_id']
        
        # handle error
        if not self.json2hdf5(source_file_id, targz_file, hdf5_file):
            output_metadata.set_exception(
                Exception(
                    "json2hdf5: Could not process files {}, {}.".format(*input_files)))
        
        return ([hdf5_file], [output_metadata])
