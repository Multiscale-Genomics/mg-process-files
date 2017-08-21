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

from __future__ import print_function

import argparse
import os

# Required for ReadTheDocs
from functools import wraps # pylint: disable=unused-import

from basic_modules.workflow import Workflow
from dmp import dmp

from tool.json_3d_indexer import json3dIndexerTool

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


    def run(self, input_files, metadata, output_files):
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

        json_tar_file = input_files[0]
        hdf5_file = input_files[1]

        # 3D JSON Indexer
        j = json3dIndexerTool()
        hdf5_idx, hdf5_meta = j.run([json_tar_file, hdf5_file], [], metadata)

        return (hdf5_idx, hdf5_meta)

# ------------------------------------------------------------------------------

def main(input_files, output_files, input_metadata):
    """
    Main function
    -------------

    This function launches the app.
    """

    # import pprint  # Pretty print - module for dictionary fancy printing

    # 1. Instantiate and launch the App
    print("1. Instantiate and launch the App")
    from apps.workflowapp import WorkflowApp
    app = WorkflowApp()
    result = app.launch(process_json_3d, input_files, input_metadata,
                        output_files, {})

    # 2. The App has finished
    print("2. Execution finished")
    print(result)
    return result

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Index the bed file")
    PARSER.add_argument("--assembly", help="Assembly")
    PARSER.add_argument("--gz_file", help="Compressed tar file of 3D JSON files to get indexed")
    PARSER.add_argument("--h5_file", help="Location of HDF5 index file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    GZ_FILE = os.path.abspath(ARGS.gz_file)
    ASSEMBLY = ARGS.assembly
    HDF5_FILE = ARGS.h5_file

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    are ones that are present in the ENA then I would need to download them


    #2. Register the data with the DMP
    DM_HANDLER = dmp(test=True)

    print(DM_HANDLER.get_files_by_user("test"))

    J_FILE = DM_HANDLER.set_file(
        "test", GZ_FILE, "gz", "Model", "", "gz", [],
        {"assembly" : ASSEMBLY}
    )
    H5_FILE = DM_HANDLER.set_file(
        "test", HDF5_FILE, "hdf5", "index", "", None, [J_FILE],
        {"assembly" : ASSEMBLY}
    )

    print(DM_HANDLER.get_files_by_user("test"))

    # 3. Instantiate and launch the App
    #from basic_modules import WorkflowApp
    #app = WorkflowApp()
    #results = app.launch([J_FILE, h5_file], {"assembly" : assembly})

    RESULTS = main(
        [GZ_FILE, HDF5_FILE], [],
        {"assembly" : ASSEMBLY, "file_id" : J_FILE}
    )

    print(DM_HANDLER.get_files_by_user("test"))
