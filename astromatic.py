#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 14:05:34 2024

@author: amartinez
"""
import subprocess
import sys
import time
import glob
import os
import shutil
# This scrip run all the three astromatic programs used for the corrections of geometric distorsion in GNS1
# %%           
field = 'B6'
band = 'J'
sex_folder =   '/home/data/alvaro/gns_gd/gns1/%s/%s/sextractor/'%(band,field)
scamp_folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/scamp/'%(band,field)
SWarp_folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/SWarp/'%(band,field)
cubes_aligned ='/home/data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/'%(band,field)
scripts = '/home/data/alvaro/gns_gd/gns1/scripts/'
out_folder = SWarp_folder + 'outputs/'

# %%
ch_range = [1,5]
#SOURCE-EXTRACTOR
t0_sex = time.time()
for chip in range(ch_range[0],ch_range[1]):
    command = ['source-extractor', cubes_aligned + '%s_image_c%s.fits'%(field,chip), 
               '-c', scripts + 'default_c.sex', '-CATALOG_NAME',f'{sex_folder}chip{chip}/%s_image_c%s.cat'%(field, chip),
               '-CHECKIMAGE_NAME',f'{sex_folder}chip{chip}/objects_%s_c%s.fits'%(field, chip)]
    
    try:
        # Run the command
        
        # result = subprocess.run(command, cwd=f'{sex_folder}chip{chip}/',check=True, text=True, capture_output=True)
        # result = subprocess.run(command, cwd=f'{sex_folder}chip{chip}/',check=True)
        result = subprocess.run(command, cwd= scripts,check=True)
        
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
                '-c', 'scamp_c.conf', '-HEADER_NAME', f'{scamp_folder}chip{chip}/%s_image_c%s.head'%(field, chip)
                ,'-FULLOUTCAT_NAME',f'{scamp_folder}chip{chip}/%s_full_c%s.cat'%(field, chip)]
    
    try:
        # Run the command
        
        result = subprocess.run(command, cwd=scripts,check=True)
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
    fits_files = glob.glob(os.path.join(scripts, '*.pdf'))  # Find all .pdf files
    destination_folder = f'{scamp_folder}chip{chip}/'
    
    for fits_file in fits_files:
        try:
            # Construct destination file path
            dest_file = os.path.join(destination_folder, os.path.basename(fits_file))
            
            # Remove the file if it already exists in the destination folder
            if os.path.exists(dest_file):
                os.remove(dest_file)
                print(f"Existing file removed: {dest_file}")
            
            # Move the file
            shutil.move(fits_file, destination_folder)
            print(f"Moved: {fits_file} -> {dest_file}")
        except Exception as e:
            print(f"Error moving {fits_file}: {e}")
t1_sca = time.time()
t_sca = t1_sca-t0_sca
print(30*'_' + '\nDone with Scamp\nIt tooks %.0f sec\n'%(t_sca) + 30*'_')

# %%
#SWARP
t0_swa = time.time()
for chip in range(ch_range[0],ch_range[1]):
    
    command = ['SWarp', cubes_aligned+ '%s_image_c%s.fits' %(field, chip), 
               '-c', 'default_c.swarp', '-HEADER_NAME',scamp_folder + 'chip%s/%s_image_c%s.head'%(chip,field, chip),
               '-WEIGHT_IMAGE',cubes_aligned + '%s_mask_c%s.fits'%(field, chip),
               '-RESAMPLE_DIR', out_folder + 'chip%s/'%(chip)]
    
    try:
        # Run the command
        
        # result = subprocess.run(command, cwd=f'{SWarp_folder}chip{chip}/',check=True)
        result = subprocess.run(command, cwd= scripts ,check=True)
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
