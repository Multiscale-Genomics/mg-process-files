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

from tool.gff3_sorter import gff3SortTool
from tool.gff3_indexer import gff3IndexerTool


@pytest.mark.gff3
def test_gff3_01_sorter():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "gff3": resource_path + "sample.gff3"
    }

    output_files = {
        "sorted_gff3": resource_path + "sample.sorted.gff3"
    }

    metadata = {
        "gff3": Metadata(
            "data_rnaseq", "gff3", [], None,
            {'assembly': 'test'}),
    }

    bs_handle = gff3SortTool()
    bs_handle.run(input_files, metadata, output_files)

    print(resource_path)
    assert os.path.isfile(resource_path + "sample.sorted.gff3") is True
    assert os.path.getsize(resource_path + "sample.sorted.gff3") > 0


@pytest.mark.gff3
def test_gff3_02_indexer():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    f_check = h5py.File(resource_path + "file_index.hdf5", "a")
    f_check.close()

    input_files = {
        "gff3": resource_path + "sample.sorted.gff3",
        "chrom_file": resource_path + "chrom_GRCh38.size",
        "hdf5_file": resource_path + "file_index.hdf5"
    }

    output_files = {
        "gz_file": resource_path + "sample.gff3.gz",
        "tbi_file": resource_path + "sample.gff3.gz.tbi"
    }

    metadata = {
        "gff3": Metadata(
            "data_rnaseq", "gff3", "test_gff3_location", [], {'assembly': 'test'}),
        "hdf5_file": Metadata(
            "data_file", "hdf5", "test_location", [], {}
        )
    }

    gi_handle = gff3IndexerTool()
    gi_handle.run(input_files, metadata, output_files)

    print(resource_path)
    assert os.path.isfile(resource_path + "sample.gff3.gz") is True
    assert os.path.getsize(resource_path + "sample.gff3.gz") > 0
    assert os.path.isfile(resource_path + "sample.gff3.gz.tbi") is True
    assert os.path.getsize(resource_path + "sample.gff3.gz.tbi") > 0
