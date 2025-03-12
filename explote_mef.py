#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:18:34 2024

@author: amartinez
"""

import subprocess

# %
# %%
# %%
# #MISSFITS

files = ['clean_FAST-Bulge-6.Ks.2016-06-14T00_43_38']
command = ['missfits', '*.fits', '-SAVE_TYPE', 'BACKUP', '-OUTFILE_TYPE', 'DIR', ]
    
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





        
