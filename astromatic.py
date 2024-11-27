#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:05:34 2024

@author: amartinez
"""
import subprocess
import sys
import time
# This scrip run all the three astromatic programs used for the corrections of geometric distorsion in GNS1
# %%           
field = 'B6'
band = 'Ks'
sex_folder =   '/home/data/alvaro/gns_gd/gns1/%s/%s/sextractor/'%(band,field)
scamp_folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/scamp/'%(band,field)
SWarp_folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/SWarp/'%(band,field)
cubes_aligned ='/home/data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/'%(band,field)
out_folder = SWarp_folder + 'outputs/'
sys.exit(20)
# %%
ch_range = [1,2]
#SOURCE-EXTRACTOR
t0_sex = time.time()
for chip in range(ch_range[0],ch_range[1]):
    command = ['source-extractor', cubes_aligned + '%s_image_c%s.fits'%(field,chip), 
               '-c', 'default_c%s.sex'%(chip), '-CATALOG_NAME','%s_image_c%s.cat'%(field, chip)
               ]
    
    try:
        # Run the command
        
        # result = subprocess.run(command, cwd=f'{sex_folder}chip{chip}/',check=True, text=True, capture_output=True)
        result = subprocess.run(command, cwd=f'{sex_folder}chip{chip}/',check=True)
        
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
t1_sex = time.time()

t_sex = t1_sex-t0_sex
print(30*'_' + '\nDone with SExtractor\nIt tooks %.0f sec\n'%(t_sex) + 30*'_')
# %%
#SCAMP
t0_sca = time.time()
for chip in range(ch_range[0],ch_range[1]):
    command = ['scamp', sex_folder + 'chip%s/%s_image_c%s.cat'%(chip, field,chip), 
                '-c', 'scamp_c%s.conf'%(chip), '-HEADER_NAME', '%s_image_c%s.head'%(field, chip)
                ,'-FULLOUTCAT_NAME','%s_full_c%s.cat'%(field, chip)]
    
    try:
        # Run the command
        
        result = subprocess.run(command, cwd=f'{scamp_folder}chip{chip}/',check=True)
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
t1_sca = time.time()
t_sca = t1_sca-t0_sca
print(30*'_' + '\nDone with Scamp\nIt tooks %.0f sec\n'%(t_sca) + 30*'_')

#%%
#SWARP
t0_swa = time.time()
for chip in range(ch_range[0],ch_range[1]):
    
    command = ['SWarp', cubes_aligned+ '%s_image_c%s.fits' %(field, chip), 
               '-c', 'default_c%s.swarp'%(chip), '-HEADER_NAME',scamp_folder + 'chip%s/%s_image_c%s.head'%(chip,field, chip),
               '-WEIGHT_IMAGE',cubes_aligned + '%s_mask_c%s.fits'%(field, chip),
               '-RESAMPLE_DIR', out_folder + 'chip%s/'%(chip)]
    
    try:
        # Run the command
        
        result = subprocess.run(command, cwd=f'{SWarp_folder}chip{chip}/',check=True)
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
t1_swa = time.time()
t_swa = t1_swa-t0_swa
print(30*'_' + '\nDone with SWarp\nIt tooks %.0f sec\n'%(t_swa) + 30*'_')
