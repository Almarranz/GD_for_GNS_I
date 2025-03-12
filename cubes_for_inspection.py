#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:18:34 2024

@author: amartinez
"""

import astropy.io.fits as fits
import numpy as np
import os
import sys
from astropy.table import Table
from astropy.wcs.utils import fit_wcs_from_points
from astropy.wcs import WCS
from astropy.io import fits
import astroalign as aa
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs.utils import fit_wcs_from_points
import shutil
import gzip
import subprocess

# %
# %%
field = 11

ch_range = [2,3]


# %%
# #MISSFITS

for chip in range(ch_range[0],ch_range[1]):
    pruebas = '/home/data/alvaro/gns_gd/gns1/pruebas/part1_B1_B6/chip%s/'%(chip)
    command = [pruebas, 'missfits', pruebas + 'clean_FAST-Bulge-6.Ks.2016-06-14T00_43_38_%s'%(chip), '-SPLIT_SUFFIX', '.%03d.fits', '-SAVE_TYPE', 'REPLACE', '-OUTFILE_TYPE', 'MULTI']
    
    try:
        # Run the command
        
        result = subprocess.run(command, check=True)
        # Print standard output and error
        print("Command Output:")
        print(result.stdout)
        print("Command Error (if any):")
        print(result.stderr)
    
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Error: {e}")
        print(f"Standard Output: {e.stdout}")
        print(f"Standard Error: {e.stderr}")





        
