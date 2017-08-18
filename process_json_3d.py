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

from tool import json_3d_indexer

# ------------------------------------------------------------------------------

class process_json_3d(Workflow):
    """
    Workflow to index JSON formatted files within the Multiscale Genomics (MuG)
    Virtural Research Environment (VRE) that have been generated as part of the
    Hi-C analysis pipeline to model the 3D structure of the genome within the
    nucleus of the cell.
    """

    def __init__(self):
        """
        Initialise the class
        """


    def run(self, file_ids, metadata):
        """
        Main run function to index the 3D JSON files that have been generated
        as part of the Hi-C analysis pipeline to model the 3D structure of the
        genome within the nucleus of the cellready for use in the RESTful API.

        Parameters
        ----------
        files_ids : list
            file : str
                Location of the tar.gz file of JSON files representing the 3D
                models of the nucleus
        metadata : list

        Returns
        -------
        outputfiles : list
            List with the location of the HDF5 index file for the given dataset
        """

        json_tar_file = file_ids[0]
        h5_file = file_ids[1]

        # 3D JSON Indexer
        j = json_3d_indexer.json3dIndexerTool()
        h5_idx = j.run([json_tar_file, h5_file], metadata)

        return (h5_idx)

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    import os

    # Set up the command line parameters
    parser = argparse.ArgumentParser(description="Index the bed file")
    parser.add_argument("--assembly", help="Assembly")
    parser.add_argument("--gz_file", help="Compressed tar file of 3D JSON files to get indexed")
    parser.add_argument("--h5_file", help="Location of HDF5 index file")

    # Get the matching parameters from the command line
    args = parser.parse_args()

    gz_file = os.path.abspath(args.gz_file)
    assembly = args.assembly
    hdf5_file = os.path.abspath(args.h5_file)

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    are ones that are present in the ENA then I would need to download them


    #2. Register the data with the DMP
    da = dmp(test=True)

    print da.get_files_by_user("test")

    j_file = da.set_file("test", gz_file, "gz", "Model", "", "gz", [], {"assembly" : assembly})
    h5_file = da.set_file("test", hdf5_file, "hdf5", "index", "", None, [j_file], {"assembly" : assembly})

    print da.get_files_by_user("test")

    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch([j_file, h5_file], {"assembly" : assembly})

    pj = process_json_3d()
    results = pj.run([gz_file, hdf5_file], {"assembly" : assembly, "file_id" : j_file})

    print da.get_files_by_user("test")
