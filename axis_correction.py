# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Thu Oct 24 11:35:06 2024
@author: amartinez
"""

# %%
# Modify the size of the axes from the output of SWarp to make them equal
import numpy as np
from astropy.io import fits
import os
from astropy.table import Table
import sys
import time
import glob
pruebas = '/home/data/alvaro/gns_test/gns1/B6/geometric_dis/pruebas/'

field = 'B6'
band = 'Ks'
ax1_size = 2048
ax2_size = 768
ch_range = [1,2]
# sys.exit(24)
# %% /Volumes/teabag-data
for chip in range(ch_range[0],ch_range[1]):
    tc0 = time.time()
    #This run the program from teatime
    # folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/SWarp/outputs/chip%s/'%(band,field, chip)
    # folder_gd = '/home/data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/chip%s/'%(band, field, chip)
    # slice_list =  '/home/data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/'%(band, field)
    
    # This run the program from my local machine
    folder = '/Volumes/teabag-data/alvaro/gns_gd/gns1/%s/%s/SWarp/outputs/chip%s/'%(band,field, chip)
    folder_gd = '/Volumes/teabag-data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/chip%s/'%(band, field, chip)
    slice_list =  '/Volumes/teabag-data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/'%(band, field)
    
    # folder = pruebas
    fits_files = [f for f in sorted(os.listdir(folder)) if f.endswith('.fits') and f.startswith('%s_image_c%s' % (field, chip)) or f.startswith('%s_image_c%s' % (field, chip)) and f.endswith('resamp.weight.fits')]
    # fits_files = fits_files[0:2] 
    # sys.exit(25)    
    for nf, f_file in enumerate(fits_files):
        hdul = fits.open(folder + f_file)
        image_data = hdul[0].data
        ejes = [hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2']]
        print(f'Working with {f_file}')
        print(f'Original axes: {ejes}')
        
        # arget_size = 2048
        axis1, axis2 = ejes

        # Calculate crop/padding amounts for each axis
        crop_pad_axis1 = (axis1 - ax1_size) // 2
        crop_pad_axis2 = (axis2 - ax2_size) // 2
        
        if (axis1 - ax1_size) < -0 or (axis2 - ax2_size) < -0:
            print(30*'∫' + f'\nCheck this one {f_file}\n' + 30*'∫')

        
        # Process for axis 1 if not equal to target size
        if axis1 != ax1_size:
            if axis1 > ax1_size:  # Crop
                start_x = crop_pad_axis1
                end_x = start_x + ax1_size
                image_data = image_data[:, start_x:end_x]
                hdul[0].header['CRPIX1'] -= start_x
            else:  # Pad
                if (axis1 - ax1_size) % 2 == 0: 
                    pad_x = abs(crop_pad_axis1)
                    image_data = np.pad(image_data, ((0, 0), (pad_x, pad_x)), mode='constant', constant_values=0)
                    hdul[0].header['CRPIX1'] += pad_x
                else:
                    resto = (axis1 - ax1_size) %2
                    pad_x = abs(crop_pad_axis1)
                    image_data = np.pad(image_data, ((0, 0), (pad_x, pad_x + resto)), mode='constant', constant_values=0)
                    hdul[0].header['CRPIX1'] += pad_x

        # Process for axis 2 if not equal to target size
        if axis2 != ax2_size:
            if axis2 > ax2_size:  # Crop
                start_y = crop_pad_axis2
                end_y = start_y + ax2_size
                image_data = image_data[start_y:end_y, :]
                hdul[0].header['CRPIX2'] -= start_y
            else:  # Pad
                if (axis2 - ax2_size) %2 == 0:
                    pad_y = abs(crop_pad_axis2)
                    image_data = np.pad(image_data, ((pad_y, pad_y), (0, 0)), mode='constant', constant_values=0)
                    hdul[0].header['CRPIX2'] += pad_y
                else:
                    resto = (axis2 - ax2_size) %2
                    print('YOMAMMAAAAAAAAA')
                    pad_y = abs(crop_pad_axis2)
                    image_data = np.pad(image_data, ((pad_y, pad_y + resto), (0, 0)), mode='constant', constant_values=0)
                    hdul[0].header['CRPIX2'] += pad_y

        # Update header axes and save to file
        hdul[0].header['NAXIS1'] = ax1_size
        hdul[0].header['NAXIS2'] = ax2_size
        # hdul[0].data = image_data
        if f_file.endswith('weight.fits'):
            image_data = np.where(image_data ==0,0,1)
            hdul[0].data = image_data
        else:
            hdul[0].data = image_data
        hdul.writeto(folder + f'PADDED_{f_file}', overwrite=True)
        
        print(f'New axes: {hdul[0].header["NAXIS1"]}, {hdul[0].header["NAXIS2"]}')
        print(30 * '_')
    tc1 = time.time()   
    tc = tc1 - tc0
    print(30*'*' + f'\nDone with chip {chip}\nIt tooks {tc:.0f} sec\n' + 30*'*')

# %%

sl = Table.read(slice_list + '%s_cubes_and_slices.txt'%(field), format = 'ascii')

ax1_sz = ax1_size
ax2_sz = ax2_size

# chip = 1
for chip in range(ch_range[0],ch_range[1]):
    folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/SWarp/outputs/chip%s/'%(band,field, chip)
    folder_gd = '/home/data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/chip%s/'%(band, field, chip)
    # for n,fi in enumerate(fits_files):
    idex = 1
    for i in range(len(sl)):
        slices = sl['number_of_slices'][i]
        cube_name = sl['Cube_id'][i]
        print(sl['Cube_id'][i], sl['number_of_slices'][i])
        cube = np.empty((slices, ax2_sz, ax1_sz))
        cube_w = np.empty((slices, ax2_sz, ax1_sz))
        for j in range(slices):
            hdul = fits.open(folder + 'PADDED_%s_image_c%s.%04d.resamp.fits'%(field,chip,idex))
            hdul_w = fits.open(folder + 'PADDED_%s_image_c%s.%04d.resamp.weight.fits'%(field,chip,idex))
            print('PADDED_%s_image_c%s.%04d.resamp.fits'%(field,chip,idex))
            cube[j,:,:] = hdul[0].data    
            cube_w[j,:,:] = hdul_w[0].data    
            idex +=1
        fits.writeto(folder_gd + 'cube%s_gd.fits'%(cube_name),cube, overwrite=True)
        fits.writeto(folder_gd + 'cube%s_gd_w.fits'%(cube_name),cube_w, overwrite=True)
        
# # %
# =============================================================================
# # This part delete de padded and intermediate products
# =============================================================================        
    
#       # Iterate over each chip subfolder
# for chip in range(ch_range[0],ch_range[1]):
#     folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/SWarp/outputs/chip%s/'%(band,field, chip)
#     fits_files = glob.glob(os.path.join(folder, '*.fits'))  # Find all .fits files
#     # for i in fits_files:
#     #     print(i)
#     for fits_file in fits_files:
#         try:
#             os.remove(fits_file)  # Remove the FITS file
#             print(f"Deleted: {fits_file}")
#         except Exception as e:
#             print(f"Error deleting {fits_file}: {e}")

# print("All FITS files is SWarp removed.")  
        
# for chip in range(ch_range[0],ch_range[1]):
#     folder = '/home/data/alvaro/gns_gd/gns1/%s/%s/sextractor/chip%s/'%(band,field, chip)
#     fits_files = glob.glob(os.path.join(folder, '*.fits'))  # Find all .fits files
#     # for i in fits_files:
#     #     print(i)
#     for fits_file in fits_files:
#         try:
#             os.remove(fits_file)  # Remove the FITS file
#             print(f"Deleted: {fits_file}")
#         except Exception as e:
#             print(f"Error deleting {fits_file}: {e}")
        
        
# print("All FITS files is sextractor removed.")  
# %%
# =============================================================================
# # Be careful!!This part delete the original aligned cubes.
# =============================================================================

# =============================================================================
#   # Assuming subfolders are named chip1, chip2, ..., chip4
# folder =  '/home/data/alvaro/gns_gd/gns1/%s/%s/cubes_gd/'%(band, field)
# fits_files = glob.glob(os.path.join(folder, '*.fits'))  # Find all .fits files
# for i in fits_files:
#     print(i)
# for fits_file in fits_files:
#     try:
#         os.remove(fits_file)  # Remove the FITS file
#         print(f"Deleted: {fits_file}")
#     except Exception as e:
#         print(f"Error deleting {fits_file}: {e}")
#          
#         
# print("All FITS aligned cubes removed.")  
# =============================================================================
# %%




