#!/bin/python

import os
from glob import glob
from bruker2nifti.converter import Bruker2Nifti

data_in = '/home/julia/projects/real_data/mouse_visual/CL181102fmrssouris9/raw/'
data_out = '/home/julia/projects/real_data/mouse_visual/CL181102fmrssouris9/raw/26/converted/'

if not os.path.isdir(data_out):
    os.mkdir(data_out)

bru = Bruker2Nifti(data_in, data_out)
bru.get_acqp = True
bru.get_method = True
bru.get_reco = True
bru.verbose = 0
bru.convert()
