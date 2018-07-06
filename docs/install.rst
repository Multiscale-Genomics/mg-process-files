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

Requirements and Installation
=============================

Requirements
------------

Software
^^^^^^^^

- Mongo DB 3.2
- Python 2.7.10+
- SamTools
- bedToBigBed - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- wigToBigWig - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- HDF5

Python Modules
^^^^^^^^^^^^^^

- numpy
- h5py
- pyBigWig
- pysam

Installation
------------
Directly from GitHub:

.. code-block:: none
   :linenos:

   git clone https://github.com/Multiscale-Genomics/mg-process-files.git

Using pip:

.. code-block:: none
   :linenos:

   pip install git+https://github.com/Multiscale-Genomics/mg-process-files.git