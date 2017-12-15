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

import sys
import subprocess
import shlex

from utils import logger

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("          Using mock decorators.")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task
    from utils.dummy_pycompss import compss_wait_on

from basic_modules.metadata import Metadata
from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class bedSortTool(Tool):
    """
    Tool for running indexers over a BED file for use in the RESTful API
    """

    def __init__(self, configuration=None):
        """
        Init function
        """
        logger.info("BED File Sorter")
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)

    @task(returns=bool, bed_file=FILE_IN, bed_out_file=FILE_OUT)
    def bedsorter(self, file_bed, bed_out_file): # pylist disable=could-be-function
        """
        BED file sorter

        This is a wrapper for the standard Linux ``sort`` method the sorting by
        the chromosome and start columns in the BED file.

        Parameters
        ----------
        file_bed : str
            Location of the BED file to sort

        Example
        -------
        .. code-block:: python
           :linenos:

           results = self.bedsorter(bed_file)
           results = compss_wait_on(results)

        """
        tmp_bed_file = file_bed + '.sorted.bed'
        with open(tmp_bed_file, "wb") as out:
            command_line = 'sort -k1,1 -k2,2n ' + file_bed
            args = shlex.split(command_line)
            process_handle = subprocess.Popen(args, stdout=out)
            process_handle.wait()

        with open(bed_out_file, "wb") as f_out:
            with open(tmp_bed_file, "rb") as f_in:
                f_out.write(f_in.read())

        return True


    def run(self, input_files, input_metadata, output_files):
        """
        Function to run the BED file sorter

        Parameters
        ----------
        input_files : list
            bed_file : str
                Location of the bed file
        metadata : list
            file_id : str
                file_id used to identify the original bed file
            assembly : str
                Genome assembly accession

        Returns
        -------
        list
            bed_file : str
                Location of the sorted bed file

        Example
        -------
        .. code-block:: python
           :linenos:

           from tool.bed_sorter import bedSortTool

           # Bed Sorter
           bst = bedSortTool()
           bst_files, bst_meta = bst.run([bed_file], [], {})
        """
        results = self.bedsorter(input_files["bed"], output_files["sorted_bed"])
        results = compss_wait_on(results)

        output_metadata = {
            "sorted_bed": Metadata(
                data_type=input_metadata["bed"].data_type,
                file_type=input_metadata["bed"].file_type,
                file_path=input_metadata["bed"].file_path,
                sources=[],
                taxon_id=input_metadata["bed"].taxon_id,
                meta_data={
                    "tool": "bed_sorter"
                }
            )
        }

        return ({"sorted_bed" : input_files["bed"]}, output_metadata)

# ------------------------------------------------------------------------------
