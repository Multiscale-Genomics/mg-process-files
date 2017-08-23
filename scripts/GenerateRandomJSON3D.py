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

import numpy as np
import json

meta_data = {
    "version"  : 1.0,
    "type"     : "dataset",
    "generator": "TADbit"
}

objects = {
    "chromEnd" : [20440000.0],
    "chromStart" : [20420000.0],
    "end" : 22000,
    "chrom" : ["chr2R"],
    "experimentType" : "Hi-C",
    "species" : "Drosophila melanogaster",
    "project" : "Test_dataset",
    "start" : 2000,
    "assembly" : "DBGP 5 (dm3)",
    "identifier" : "SRRnnnnn",
    "resolution" : 2000,
    "cellType" : "kc167",
    "uuid": "unique_id_1234567890",
    "title": "Sample TADbit data",
    "datatype": "xyz",
    "components": 3,
    "source": "local",
    "dependencies": {
        "IMP": "2.0.1 (random seed indexed at 1 = 89400484)",
        "TADbit": "0.1_alpha.723",
        "MCL": "12-135"
    }
}

models = []
hic_data = {"data" : {}, "n" : 21, "tads" : [[[1,20420000,20430000,10], [1,20430000,20440000,10]]]}

for i in range(10):
    data = list(np.random.randint(20, size=21)-10)
    models.append({"ref" : i, "data" : data})

r_data = list(np.random.random_sample((440,)))
for i in range(440):
    hic_data["data"][str(i)] = r_data[i]

clusters = [[1,3,5,7,9], [2,4,6,8,0]]
centroids = [1, 6]

restraints = []

json_models = {
    "metadata" : meta_data,
    "object" : objects,
    "models" : models,
    "clusters" : clusters,
    "centroids" : centroids,
    "restraints" : restraints,
    "hic_data" : hic_data
}

f = open("../tests/data/sample_3D_models.json", "w")
f.write(json.dumps(json_models))
f.close()