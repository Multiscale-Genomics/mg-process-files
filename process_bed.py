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
import h5py

from basic_modules.workflow import Workflow
from utils import logger
from utils import remap

from tool.bed_indexer import bedIndexerTool
from tool.bed_sorter import bedSortTool

# ------------------------------------------------------------------------------

class process_bed(Workflow):
    """
    Workflow to index BED formatted files within the Multiscale Genomics (MuG)
    Virtural Research Environment (VRE)
    """

    def __init__(self, configuration=None):
        """
        Initialise the tool with its configuration.


        Parameters
        ----------
        configuration : dict
            a dictionary containing parameters that define how the operation
            should be carried out, which are specific to each Tool.
        """
        logger.info("Process BED file")
        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    def run(self, input_files, metadata, output_files):
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
        inpout_files : list
            List of file locations
        metadata : list

        Returns
        -------
        outputfiles : list
            List of locations for the output BED and HDF5 files
        """

        # Ensure that the file exists
        f_check = h5py.File(input_files["hdf5_file"], "a")
        f_check.close()

        # Bed Sort
        bst = bedSortTool()
        bst_files, bst_meta = bst.run(
            remap(input_files,
                  "bed"),
            remap(metadata,
                  "bed"),
            remap(output_files,
                  "sorted_bed")
        )

        # Bed Indexer
        bit = bedIndexerTool()
        bit_files, bit_meta = bit.run(
            {
                "bed" : bst_files["sorted_bed"],
                "chrom_size" : input_files["chrom_size"],
                "hdf5_file" : input_files["hdf5_file"]
            }, {
                "bed" : bst_meta["sorted_bed"],
                "chrom_size" : metadata["chrom_size"],
                "hdf5_file" : metadata["hdf5_file"]
            }, {
                "bb_file" : output_files["bb_file"],
                "hdf5_file" : output_files["hdf5_file"]
            }
        )

        return (bst_files + bit_files, bst_meta + bit_meta)

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
    result = app.launch(process_bed,
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
    print(RESULTS)
