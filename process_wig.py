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

import h5py

from basic_modules.workflow import Workflow
from dmp import dmp

from tool.wig_indexer import wigIndexerTool


# ------------------------------------------------------------------------------

class process_wig(Workflow):
    """
    Workflow to index WIG formatted files within the Multiscale Genomics (MuG)
    Virtural Research Environment (VRE)
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
        if configuration is None:
            configuration = {}
        self.configuration.update(configuration)


    def run(self, input_files, metadata, output_files):
        """
        Main run function to index the WIG files ready for use in the RESTful
        API. WIG files are indexed in 2 different ways to allow for optimal data
        retreival. The first is as a bigwig file, this allows the data to get
        easily extracted as WIG documents and served to the user. The second is
        as an HDF5 file that is used to identify which bed files have
        information at a given location. This is to help the REST clients make
        only the required calls to the relevant WIG files rather than needing to
        pole all potential WIG files.

        Parameters
        ----------
        files_ids : list
            List of file locations
        metadata : list

        Returns
        -------
        outputfiles : list
            List of locations for the output wig and HDF5 files
        """

        wig_file = input_files[0]
        chrom_file = input_files[1]
        hdf5_file = input_files[2]

        # Ensure that the file exists
        f_check = h5py.File(hdf5_file, "a")
        f_check.close()

        # Bed Indexer
        wit = wigIndexerTool()
        wit_files, wit_meta = wit.run((wig_file, chrom_file, hdf5_file), metadata)

        return (wit_files, wit_meta)

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
    result = app.launch(process_wig, input_files, input_metadata,
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
    PARSER.add_argument("--chrom", help="Matching chrom.size file")
    PARSER.add_argument("--wig_file", help="WIG file to get indexed")
    PARSER.add_argument("--h5_file", help="HDF5 index file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    ASSEMBLY = ARGS.assembly
    CHROM_SIZE_FILE = ARGS.chrom
    WIG_FILE = ARGS.wig_file
    HDF5_FILE = ARGS.h5_file

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    Are ones that are present in the ENA then I would need to download them


    #2. Register the data with the DMP
    DM_HANDLER = dmp(test=True)

    print(DM_HANDLER.get_files_by_user("test"))

    CS_FILE = DM_HANDLER.set_file(
        "test", CHROM_SIZE_FILE, "tsv", "ChIP-seq", "", None, [],
        {"assembly" : ASSEMBLY}
    )
    W_FILE = DM_HANDLER.set_file(
        "test", WIG_FILE, "wig", "Assembly", "", None, [],
        {"assembly" : ASSEMBLY}
    )
    H5_FILE = DM_HANDLER.set_file(
        "test", HDF5_FILE, "hdf5", "index", "", None, [CS_FILE, W_FILE],
        {"assembly" : ASSEMBLY}
    )

    print(DM_HANDLER.get_files_by_user("test"))

    METADATA = {
        "file_id" : W_FILE,
        "assembly" : ASSEMBLY
    }

    # 3. Instantiate and launch the App
    RESULTS = main([WIG_FILE, CHROM_SIZE_FILE, HDF5_FILE], [], METADATA)

    print(DM_HANDLER.get_files_by_user("test"))
