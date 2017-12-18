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
import pytest # pylint: disable=unused-import

from basic_modules.metadata import Metadata

from tool.json_3d_indexer import json3dIndexerTool

@pytest.mark.json3d
def test_json3d_indexer():
    """
    Function to test Kallisto indexer
    """
    resource_path = os.path.join(os.path.dirname(__file__), "data/")

    input_files = {
        "models" : resource_path + "sample_3D_models.tar.gz"
    }

    output_files = {
        "index" : resource_path + "sample.models.hdf5"
    }

    metadata = {
        "models" : Metadata(
            "data_rnaseq", "gff3", "test_gff3_location", [], {'assembly' : 'test'})
    }

    j3d_handle = json3dIndexerTool()
    j3d_handle.run(input_files, metadata, output_files)

    print(resource_path)
    # assert os.path.isfile(resource_path + "sample.gff3.gz") is True
    # assert os.path.getsize(resource_path + "sample.gff3.gz") > 0
