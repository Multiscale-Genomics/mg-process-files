# mg-process-tsv

[![Documentation Status](https://readthedocs.org/projects/mg-process-tsv/badge/?version=latest)](http://mg-process-files.readthedocs.org/en/latest/)

Scripts for processing text files ready for indexing for the
RESTful servers as pipelines within the MuG VRE.

Full documentation can be found on [ReadTheDocs](http://mg-process-files.readthedocs.io)

# Requirements
- bedToBigBed - http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
- HDF5
- Python 2.7.10+
- Python Modules:
  - numpy
  - h5py

# Installation
Cloneing from GitHub:
```
git clone https://github.com/Multiscale-Genomics/mg-process-tsv.git
```
To get this to be picked up by pip if part of a webserver then:
```
pip install --editable .
```
This should install the required packages listed in the `setup.py` script.


Installation via pip:
```
pip install git+https://github.com/Multiscale-Genomics/mg-process-tsv.git
```

# Configuration file
Requires a file with the name `dmp.cnf` with the following parameters to define the MongoDB server:
```
[dmp]
host = localhost
port = 27017
user = 
pass = 
db = dmp
ftp_root = ftp://ftp.<url_root>
```
