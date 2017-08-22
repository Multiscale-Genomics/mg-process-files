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

from tool.bed_indexer import bedIndexerTool
from tool.bed_sorter import bedSortTool

# ------------------------------------------------------------------------------

class process_bed(Workflow):
    """
    Workflow to index BED formatted files within the Multiscale Genomics (MuG)
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

        bed_file = input_files[0]
        chrom_file = input_files[1]
        hdf5_file = input_files[2]
        assembly = metadata["assembly"]
        bed_type = metadata["bed_type"]
        file_id = metadata["file_id"]

        meta_data = {
            "file_id": file_id,
            "bed_type": bed_type,
            "assembly": assembly
        }

        # Ensure that the file exists
        f_check = h5py.File(hdf5_file, "a")
        f_check.close()

        # Bed Sort
        bst = bedSortTool()
        bst_files, bst_meta = bst.run([bed_file], [], {})
        print('PROCESS:', bst_files)

        # Bed Indexer
        bit = bedIndexerTool()
        bit_files, bit_meta = bit.run([bst_files[0], chrom_file, hdf5_file], [], meta_data)

        return (bst_files + bit_files, bst_meta + bit_meta)

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
    result = app.launch(process_bed, input_files, input_metadata,
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
    PARSER.add_argument("--bed_file", help="Bed file to get indexed")
    PARSER.add_argument("--h5_file", help="Location of HDF5 index file")
    PARSER.add_argument("--bed_type", help="Type of Bed file bedN[+P]", default=None)

    # Get the matching parameters from the command line
    ARGS = PARSER.parse_args()

    ASSEMBLY = ARGS.assembly
    CHROM_SIZE_FILE = ARGS.chrom
    BED_FILE = ARGS.bed_file
    HDF5_FILE = ARGS.h5_file
    BED_TYPE = ARGS.bed_type



    #
    # MuG Tool Steps
    # --------------
    #
    # 1. Create data files
    #    This should have already been done by the VRE - Potentially. If these
    #    Are ones that are present in the ENA then I would need to download them


    #2. Register the data with the DMP
    DM_HANDLER = dmp(test=True)

    CS_FILE = DM_HANDLER.set_file(
        "test", CHROM_SIZE_FILE, "tsv", "ChIP-seq", "", None, [],
        {"assembly" : ASSEMBLY}
    )
    B_FILE = DM_HANDLER.set_file(
        "test", BED_FILE, "bed", "Assembly", "", None, [],
        {"assembly" : ASSEMBLY}
    )
    H5_FILE = DM_HANDLER.set_file(
        "test", HDF5_FILE, "hdf5", "index", "", None, [CS_FILE, B_FILE],
        {"assembly" : ASSEMBLY}
    )

    print(DM_HANDLER.get_files_by_user("test"))

    # 3. Instantiate and launch the App

    METADATA = {
        "file_id" : B_FILE,
        "bed_type" : BED_TYPE,
        "assembly" : ASSEMBLY
    }

    RESULTS = main([BED_FILE, CHROM_SIZE_FILE, HDF5_FILE], [], METADATA)

    print(RESULTS)
    print(DM_HANDLER.get_files_by_user("test"))
