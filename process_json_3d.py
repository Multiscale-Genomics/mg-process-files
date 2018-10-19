#!/usr/bin/python

"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

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

from basic_modules.workflow import Workflow
from utils import logger

from mg_process_files.tool.json_3d_indexer import json3dIndexerTool

# ------------------------------------------------------------------------------


class process_json_3d(Workflow):
    """
    Workflow to index JSON formatted files within the Multiscale Genomics (MuG)
    Virtural Research Environment (VRE) that have been generated as part of the
    Hi-C analysis pipeline to model the 3D structure of the genome within the
    nucleus of the cell.
    """

    configuration = {}

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("Process 3D Models file")

        if configuration is None:
            configuration = {}
        self.configuration.update(configuration)

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

        j3dit = json3dIndexerTool()
        hdf5_idx, hdf5_meta = j3dit.run(
            {
                "modles": input_files["models"],
            }, {
                "models": metadata["models"],
            }, {
                "index": output_files["index"]
            }
        )

        return (hdf5_idx, hdf5_meta)

# ------------------------------------------------------------------------------


def main_json(config, in_metadata, out_metadata):
    """
    Alternative main function
    -------------

    This function launches the app using configuration written in
    two json files: config.json and input_metadata.json.
    """
    # 1. Instantiate and launch the App
    logger.info("1. Instantiate and launch the App")
    from apps.jsonapp import JSONApp
    app = JSONApp()
    result = app.launch(process_json_3d,
                        config,
                        in_metadata,
                        out_metadata)

    # 2. The App has finished
    logger.info("2. Execution finished; see " + out_metadata)

    return result

# ------------------------------------------------------------------------------


if __name__ == "__main__":
    import sys
    sys._run_from_cmdl = True  # pylint: disable=protected-access

    # Set up the command line parameters
    PARSER = argparse.ArgumentParser(description="Index BED file")
    PARSER.add_argument("--config", help="Configuration file")
    PARSER.add_argument("--in_metadata", help="Location of input metadata file")
    PARSER.add_argument("--out_metadata", help="Location of output metadata file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    CONFIG = ARGS.config
    IN_METADATA = ARGS.in_metadata
    OUT_METADATA = ARGS.out_metadata

    RESULTS = main_json(CONFIG, IN_METADATA, OUT_METADATA)
