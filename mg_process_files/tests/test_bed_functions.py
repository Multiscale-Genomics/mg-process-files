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

import os.path
import h5py
import pytest  # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from mg_process_files.tool.bed_sorter import bedSortTool
from mg_process_files.tool.bed_indexer import bedIndexerTool


@pytest.mark.bed
def test_bed_01_sorter():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "bed": resource_path + "sample.bed"
    }

    output_files = {
        "sorted_bed": resource_path + "sample.sorted.bed"
    }

    metadata = {
        "bed": Metadata(
            "data_rnaseq", "bed", [], None,
            {'assembly': 'test'}),
    }

    bs_handle = bedSortTool()
    bs_handle.run(input_files, metadata, output_files)

    print(resource_path)
    assert os.path.isfile(resource_path + "sample.sorted.bed") is True
    assert os.path.getsize(resource_path + "sample.sorted.bed") > 0


@pytest.mark.bed
def test_bed_02_indexer():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    f_check = h5py.File(resource_path + "file_index.hdf5", "a")
    f_check.close()

    input_files = {
        "bed": resource_path + "sample.sorted.bed",
        "chrom_file": resource_path + "chrom_GRCh38.size",
        "hdf5_file": resource_path + "file_index.hdf5"
    }

    output_files = {
        "bb_file": resource_path + "sample.bb"
    }

    metadata = {
        "bed": Metadata(
            "data_rnaseq", "bed", "test_bed_location", [], {'assembly': 'test'}),
        "hdf5_file": Metadata(
            "data_file", "hdf5", "test_location", [], {}
        )
    }

    bs_handle = bedIndexerTool({"bed_type": "bed6+4"})
    bs_handle.run(input_files, metadata, output_files)

    print(resource_path)
    assert os.path.isfile(resource_path + "sample.bb") is True
    assert os.path.getsize(resource_path + "sample.bb") > 0
