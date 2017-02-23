#!/usr/bin/python

"""
Copyright 2017 EMBL-European Bioinformatics Institute

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

# chr X file X posn

import numpy as np
import h5py, os

MAX_FILES = 1024
MAX_CHROMOSOMES = 1024
MAX_CHROMOSOME_SIZE = 2000000000

user_id = 'test_user'
# Mouse
assembly = 'GCA_000001635.7'

# This is for use during testing
fi_name = 'DRR000386.sorted.bed'

resource_package = __name__
resource_path = os.path.join(os.path.dirname(__file__), 'region_idx.hdf5')

for i in range(10):
    f = h5py.File(resource_path, "a")
    file_id = 'test_file_' + '{0:03d}'.format(i)
    print('Loading', file_id)
    if str(user_id) in f:
        grp = f[str(user_id)]
        meta      = grp['meta']
        
        if str(assembly) not in grp:
            asmgrp = grp.create_group(str(assembly))
            dset = asmgrp.create_dataset('data', (0, 1, MAX_CHROMOSOME_SIZE), maxshape=(MAX_CHROMOSOMES, MAX_FILES, MAX_CHROMOSOME_SIZE), dtype='bool', chunks=True, compression="gzip")
            
            dtf = h5py.special_dtype(vlen=str)
            dtc = h5py.special_dtype(vlen=str)
            fset = asmgrp.create_dataset('files', (MAX_FILES,), dtype=dtf)
            cset = asmgrp.create_dataset('chromosomes', (MAX_CHROMOSOMES,), dtype=dtc)
            
            file_idx  = [file_id]
            chrom_idx = []
        else:
            asmgrp = grp[str(assembly)]
            dset = asmgrp['data']
            fset = asmgrp['files']
            cset = asmgrp['chromosomes']
            file_idx  = [f for f in fset if f != '']
            if file_id not in file_idx:
                file_idx.append(file_id)
                dset.resize((dset.shape[0], dset.shape[1]+1, MAX_CHROMOSOME_SIZE))
            chrom_idx = [c for c in cset if c != '']
            
    else:
        # Create the initial dataset with minimum values
        grp    = f.create_group(str(user_id))
        meta   = grp.create_group('meta')
        asmgrp = grp.create_group(str(assembly))
        
        dtf = h5py.special_dtype(vlen=str)
        dtc = h5py.special_dtype(vlen=str)
        fset = asmgrp.create_dataset('files', (MAX_FILES,), dtype=dtf)
        cset = asmgrp.create_dataset('chromosomes', (MAX_CHROMOSOMES,), dtype=dtc)
        
        file_idx  = [file_id]
        chrom_idx = []
        
        dset = asmgrp.create_dataset('data', (0, 1, MAX_CHROMOSOME_SIZE), maxshape=(MAX_CHROMOSOMES, MAX_FILES, MAX_CHROMOSOME_SIZE), dtype='bool', chunks=True, compression="gzip")
    
    # Save the list of files
    fset[0:len(file_idx)] = file_idx
    
    file_chrom_count = 0

    dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='bool')

    previous_chrom = ''
    previous_start = 0
    previous_end   = 0

    print('Row prepared ...')
    
    loaded = False
    
    with open(fi_name, 'r') as fi
        for line in fi:
            line = line.strip()
            l = line.split("\t")
            
            c = str(l[0])
            s = int(l[1])
            e = int(l[2])
            #print("Here 00")
            
            loaded = False
            
            # Need to ensure that the bed file has already been sorted
            if c != previous_chrom and previous_chrom != '':
                file_chrom_count += 1
                if previous_chrom not in chrom_idx:
                    chrom_idx.append(previous_chrom)
                    cset[0:len(chrom_idx)] = chrom_idx
                    dset.resize((dset.shape[0]+1, dset.shape[1], MAX_CHROMOSOME_SIZE))
                    #chrom_idx.append(previous_chrom)
                
                print(fi_name, str(file_idx), str(file_idx.index(file_id)), str(dset.shape), str(len(dnp)))
                dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp
                print("Loaded")
                loaded = True
                
                if file_chrom_count == 5:
                    break
                
                dnp = np.zeros([MAX_CHROMOSOME_SIZE], dtype='bool')
            #print("Here 01")
            previous_chrom = c
            dnp[s:e+1] = 1

        if loaded == False:
            if previous_chrom not in chrom_idx:
                chrom_idx.append(c)
                cset[0:len(chrom_idx)] = chrom_idx
                dset.resize((dset.shape[0]+1, dset.shape[1], MAX_CHROMOSOME_SIZE))
            
            print(fi_name, str(file_idx), str(file_idx.index(file_id)), str(dset.shape), str(len(dnp)))
            dset[chrom_idx.index(previous_chrom), file_idx.index(file_id), :] = dnp
            print("Loaded - Final B")
    
    
    fi.close()

