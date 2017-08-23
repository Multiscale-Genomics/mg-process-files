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

import sys
import subprocess
import shlex

try:
    if hasattr(sys, '_run_from_cmdl') is True:
        raise ImportError
    from pycompss.api.parameter import FILE_INOUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    print("[Warning] Cannot import \"pycompss\" API packages.")
    print("          Using mock decorators.")

    from dummy_pycompss import FILE_INOUT
    from dummy_pycompss import task
    from dummy_pycompss import compss_wait_on

from basic_modules.tool import Tool

# ------------------------------------------------------------------------------

class bedSortTool(Tool):
    """
    Tool for running indexers over a BED file for use in the RESTful API
    """

    def __init__(self):
        """
        Init function
        """
        print("BED File Sorter")
        Tool.__init__(self)


    @task(bed_file=FILE_INOUT)
    def bedsorter(self, file_bed): # pylist disable=could-be-function
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

        with open(file_bed, "wb") as f_out:
            with open(tmp_bed_file, "rb") as f_in:
                f_out.write(f_in.read())

        return True


    def run(self, input_files, output_files, metadata=None):
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
        bed_file = input_files[0]

        output_metadata = {}

        results = self.bedsorter(bed_file)
        results = compss_wait_on(results)

        return ([bed_file], [output_metadata])

# ------------------------------------------------------------------------------
