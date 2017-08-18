#(grep ^"#" in.gff3; grep -v ^"#" in.gff3 | sort -k1,1 -k4,4n) > out.sorted.gff3

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

import argparse, gzip, shutil, shlex, subprocess, os.path, json
from functools import wraps


from basic_modules import Workflow
from dmp import dmp

from tool.gff3_indexer import gff3IndexerTool
from tool.gff3_sorter import gff3SortTool

# ------------------------------------------------------------------------------

class process_gff3(Workflow):
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


    def run(self, file_ids, metadata, output_files):
        """
        Main run function to index the WIG files ready for use in the RESTful
        API. WIG files are indexed in 2 different ways to allow for optimal data
        retreival. The first is as a Tabix file, this allows the data to get
        easily extracted as GFF3 documents and served to the user. The second is
        as an HDF5 file that is used to identify which bed files have
        information at a given location. This is to help the REST clients make
        only the required calls to the relevant GFF3 files rather than needing
        to pole all potential GFF3 files.

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

        gff3_file = file_ids[0]
        chrom_file = file_ids[1]
        hdf5_file = file_ids[2]
        assembly = metadata["assembly"]

        # GFF3 Sorter
        gst = gff3SortTool()
        gst_files, gst_meta = gst.run([gff3_file], [], {})


        # GFF3 Indexer
        git = gff3IndexerTool()
        git_files, git_meta = git.run([gst_files[0], chrom_file, hdf5_file], [], {'assembly' : assembly})

        return (bb, h5_idx)

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
    result = app.launch(process_gff3, input_files, input_metadata,
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
    PARSER = argparse.ArgumentParser(description="Index the gff3 file")
    PARSER.add_argument("--assembly", help="Assembly")
    PARSER.add_argument("--gff3_file", help="GFF3 file to get indexed")
    PARSER.add_argument("--h5_file", help="HDF5 index file")

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    ASSEMBLY = ARGS.assembly
    CHROM_SIZE_FILE = ARGS.chrom
    GFF3_FILE = ARGS.gff3_file
    HDF5_FILE = ARGS.h5_file

    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    Are ones that are present in the ENA then I would need to download them


    #2. Register the data with the DMP
    DM_HANDLER = dmp()

    print(DM_HANDLER.get_files_by_user("test"))

    G3_FILE = DM_HANDLER.set_file("test", GFF3_FILE, "gff3", "Assembly", "", None)
    H5_FILE = DM_HANDLER.set_file("test", HDF5_FILE, "hdf5", "index", "", None)

    print(DM_HANDLER.get_files_by_user("test"))

    # 3. Instantiate and launch the App
    RESULTS = main([GFF3_FILE, CHROM_SIZE_FILE, HDF5_FILE], [], {"assembly" : ASSEMBLY})

    print(DM_HANDLER.get_files_by_user("test"))
