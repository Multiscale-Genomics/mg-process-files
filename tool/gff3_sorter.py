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

class gff3SortTool(Tool):
    """
    Tool for running indexers over a WIG file for use in the RESTful API
    """

    def __init__(self):
        """
        Init function
        """
        print("GFF3 File Indexer")
        Tool.__init__(self)


    @task(file_gff3=FILE_INOUT)
    def gff3Sorter(self, file_gff3):
        """
        Sorts the GFF3 file in preparation for compression and indexing

        Parameters
        ----------
        file_gff3 : str
            Location of the source GFF3 file

        GFF3 file sorter

        This is a wrapper for the standard Linux ``sort`` method the sorting by
        the chromosome and start columns in the GFF3 file.

        Parameters
        ----------
        file_gff3 : str
            Location of the source GFF3 file

        Example
        -------
        .. code-block:: python
           :linenos:

           results = self.gff3sorter(gff3_file)
           results = compss_wait_on(results)
        """
        grep_comment_command = 'grep ^"#" ' + file_gff3
        grep_gff_command = 'grep -v ^"#" ' + file_gff3
        sort_command = 'sort -k1,1 -k4,4n'

        grep_comment_args = shlex.split(grep_comment_command)
        grep_args = shlex.split(grep_gff_command)
        sort_args = shlex.split(sort_command)

        file_sorted_gff3 = file_gff3 + '.sorted.gff3'

        with open(file_sorted_gff3, "wb") as out:
            grep = subprocess.Popen(grep_comment_args, stdout=out)
            grep.wait()

        with open(file_sorted_gff3, "a") as out:
            grep = subprocess.Popen(grep_args, stdout=subprocess.PIPE)
            process_handle = subprocess.Popen(sort_args, stdin=grep.stdout, stdout=out)
            process_handle.wait()

        with open(file_gff3, "wb") as f_out:
            with open(file_sorted_gff3, "rb") as f_in:
                f_out.write(f_in.read())

        return file_gff3


    def run(self, input_files, output_files, metadata=None):
        """
        Function to run the GFF3 file sorter

        Parameters
        ----------
        input_files : list
            GFF3_file : str
                Location of the GFF3 file
        metadata : list
            file_id : str
                file_id used to identify the original GFF3 file
            assembly : str
                Genome assembly accession

        Returns
        -------
        list
            gff3_file : str
                Location of the sorted GFF3 file

        Example
        -------
        .. code-block:: python
           :linenos:

           from tool.gff3_sorter import GFF3SortTool

           # GFF3 Sorter
           bst = GFF3SortTool()
           bst_files, bst_meta = bst.run([GFF3_file], [], {})
        """
        gff3_file = input_files[0]

        output_metadata = {}

        results = self.gff3Sorter(gff3_file)
        results = compss_wait_on(results)

        return ([gff3_file], [output_metadata])

# ------------------------------------------------------------------------------
