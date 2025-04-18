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
import glob

# %
# %%
field = 'B6'


folder = '/home/data/alvaro/gns_gd/gns1/H/%s/raw/'%(field)
# cubes_aligned = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/'%(field)


pruebas = '/home/data/alvaro/gns_gd/gns1/H/%s/pruebas/'%(field)

VVV_fol = '/home/data/VVV/'
vvv_f = '/home/data/alvaro/gns_gd/gns1/VVV_fields/'
ims = '/home/data/alvaro/gns_gd/gns1/H/%s/cleaned/'%(field)
scripts = '/home/data/alvaro/gns_gd/gns1/scripts/sex_scripts/'
clean = '/home/data/alvaro/gns_gd/gns1/H/%s/cleaned/'%(field)
tmp = '/home/data/alvaro/gns_gd/gns1/H/%s/tmp/'%(field)

# field = 'B6'
# filtro = 'H'

# folder = '/home/data/raw/2015/%s/Field/%s/'%(filtro,field)
# pruebas = '/Volumes/teatime-data/alvaro/geometric_dis/%s/%s/pruebas/'%(filtro,field)
# cubes_alig = '/Volumes/teatime-data/alvaro/geometric_dis/%s/%s/cubes_aligned/'%(filtro,field)
# sf_folder = '/home/data/GNS/2015/%s/%s/data/'%(filtro,field)
# clean = '/home/data/GNS/2015/%s/%s/cleaned/'%(filtro,field)
# VVV_fol = '/home/data/working/GNS_1/VVV/'
# # ims = '/home/data/GNS/2021/H/20/ims/'
# ims = '/home/data/GNS/2015/%s/%s/ims/'%(filtro,field)
# slices_folder = '/home/data/alvaro/gns_gd/gns1/H/%s/cubes_gd/slices/'%(field)


# %%
###########IMPORTANT################
# list of raw images and raw images are in /home/data/raw/GNS_2/H/Field/20/
# list of cleand cubes is in '/home/data/GNS/2021/H/20/cleaned'
# list of starfinder stars is in '/home/data/GNS/2021/H/20/data'
# mask for all the chip is un '/home/data/GNS/2021/H/20/ims/mask.fits'
####################################

vvv_file = glob.glob(vvv_f +  f'vvv_f{field}.dat')

if len(vvv_file) > 0:
    f'Fiel {field} already covered'
else:
    vvv = Table.read(VVV_fol + 'b333-2024-05-13.dat', format = 'ascii')


lista = open(folder + 'list.txt', 'r')
lines = len(lista.readlines())# This is closing lista. I have to reopen it
lista = open(folder + 'list.txt', 'r')
# %%
hdu_m = fits.open(ims + 'mask.fits')
data_m = hdu_m[0].data
dic_sl = {}

image_i = 0
ch_range = [1,5]


last_id = 0
for li,l in enumerate(lista):
    
    print(30*'*')
    print(f'Working on file {li+1}/{lines}')
    print(30*'*')
    
    orig_header = fits.getheader(folder + l.strip())
    print(li)
    
    # In case the cubes are compressed, uncomment these lines
    # with gzip.open(clean + 'cube%s.fits.gz'%(li+1),'rb') as f_in:
    #     with open(clean + 'cube%s.fits'%(li+1),'wb') as f_out:
    #         shutil.copyfileobj(f_in, f_out)
            

    cube = fits.open(clean + 'cube%s.fits'%(li+1),mode = 'update') 
    cube[0].header = orig_header
    cube.flush()   
    image_data = cube[0].data
    
    if li>0:
        image_i = image_i + dic_sl[f'l{li-1}']
        print(30*'+',f'\nThe index is {image_i}\n',30*'+')
    
    
    
    wcs = WCS(cube[0].header, naxis = 2)
    naxis1 = cube[0].header['NAXIS1']
    naxis2 = cube[0].header['NAXIS2']
    x_off = cube[0].header['CRPIX2']
    y_off = cube[0].header['CRPIX1']
    

    print(f'Cube {li+1} has %s slices'%(cube[0].header['NAXIS3']))
    
    
    
    if len(vvv_file) < 1:
        use_idx = vvv['J']<900
        vvv_data = vvv[use_idx]
        
        vvv_ra = vvv_data['ra']
        vvv_dec = vvv_data['dec']
        # vvv_gal = SkyCoord(ra = vvv_ra, dec = vvv_dec, unit = 'degree').galactic 
        
        xy = wcs.all_world2pix(np.c_[vvv_data['ra'],vvv_data['dec']],1)
        x = xy[:,0]
        y = xy[:,1]
        
        
        vvv_data.add_columns((xy[:,0],xy[:,1]),names = ('x','y'))
        
       
        extra = 0.50
        x_new = 0
        y_new= 1500
        # mask = ((x >= -x_off) & (x < naxis1-x_off) & (y >= 0) & (y < naxis2))
        mask = ((x >= -x_off+x_new) & (x < naxis1-x_off+x_new) & (y >= y_new-int(y_new*extra)) & (y < naxis2 + y_new + int(y_new*extra)))
        
        # Apply the mask to the vvv.txt data
        vvv_overlap = vvv_data[mask]
        
        vvv_overlap.write(vvv_f  + 'vvv_f%s.dat'%(field), format = 'ascii')
    else:
        vvv_overlap = Table.read(vvv_f + f'vvv_f{field}.dat', format = 'ascii')
    # sys.exit(126)
    
    vvv_overlap.sort('J')
    
    x_vvv = vvv_overlap['x']
    y_vvv = vvv_overlap['y']
    # xh = naxis1/2
   
    
    
    xh = (max(x_vvv) + min(x_vvv))/2
    yh = (max(y_vvv) + min(y_vvv))/2
    #crop list for each chip
    for chip in range(ch_range[0],ch_range[1]):
        slices_aligned = '/home/data/alvaro/gns_gd/gns1/H/%s/cubes_gd/slices/all_chips/'%(field)

        
        if (chip == 1):
            idx = np.nonzero((x_vvv < xh) & (y_vvv < yh))
        elif (chip == 2):
            idx = np.nonzero((x_vvv > xh) & (y_vvv < yh))
        elif (chip == 3):
            idx = np.nonzero((x_vvv > xh) & (y_vvv > yh))
        elif (chip == 4):
            idx = np.nonzero((x_vvv < xh) & (y_vvv > yh))
        #
        # fig, ax = plt.subplots(1,1)
        # ax.scatter(x_vvv,y_vvv)
        # ax.scatter(x_vvv[idx],y_vvv[idx])
        # ax.axvline(xh, color = 'r')
        
        
        try:
            os.remove(tmp + f'cut_c{chip}.fits')
        except:
            print(f'NO file cut_c{chip}.fits yet')

        
        
        # try:
        #     command = ['rm', pruebas  + f'cut_c{chip}.fits']  
        #     result = subprocess.run(command, check=True)
        # except:
        #     print('NO file')
       
        if chip == 1:
            command = ['fitscopy', clean + 'cube%s.fits[1:2048,1:768,1:1]'%(li+1),tmp  + f'cut_c{chip}.fits' ]
        elif chip == 2:
            command = ['fitscopy', clean + 'cube%s.fits[2049:4096,1:768,1:1]'%(li+1),tmp  + f'cut_c{chip}.fits' ]
        elif chip == 3:
            command = ['fitscopy', clean + 'cube%s.fits[2049:4096,769:1536,1:1]'%(li+1),tmp  + f'cut_c{chip}.fits' ]
        elif chip == 4:
            command = ['fitscopy', clean + 'cube%s.fits[1:2048,769:1536,1:1]'%(li+1), tmp  + f'cut_c{chip}.fits' ]
        
        result = subprocess.run(command, check=True)
    
       
        command = ['source-extractor', tmp  + f'cut_c{chip}.fits', '-CATALOG_NAME', tmp + f'xyf_c{chip}_f{field}.cat',
                   '-PARAMETERS_NAME', scripts + 'basicASol.param']
        result = subprocess.run(command, check=True)
        os.remove(tmp + f'cut_c{chip}.fits')
        
        
        # gns = Table.read(sf_folder + 'dejitter_stars_%s_%s.txt'%(chip,li+1), format = 'ascii')
        # x_gns = gns['x']
        # y_gns = gns['y']
        
        # xy_gns = np.c_[x_gns,y_gns]
        # np.savetxt(pruebas + 'sf_f20_xy.txt',xy_gns)
        
        sources = Table.read(tmp + f'xyf_c{chip}_f{field}.cat', format = 'ascii')
        
        sources.sort('FLUX_MAX')
        sources.reverse()
        # print(sources['FLUX_MAX'][0:5])
        
        
        
        print(sources.columns)
        
        #This ara the xy coordinates from starfinder on the original cubes
        x_s = sources['X_IMAGE']
        y_s = sources['Y_IMAGE']
        
        # y_s = sources['X_IMAGE']
        # x_s = sources['Y_IMAGE']
        
        
        
        
        
        xy_s = np.c_[x_s,y_s]
        # np.savetxt(pruebas + 'sex_max_flux.txt',xy_s,fmt = '%.8f')
        
        xy_vvv = np.vstack((vvv_overlap['x'][idx],vvv_overlap['y'][idx])).T
        # Astroalign does not work when repeated value appears in the same array
        xy_test, idx_u, = np.unique(xy_vvv, axis=0, return_index=True)
        xy_vvv_unique = xy_vvv[np.sort(idx_u)]
        
        np.savetxt(pruebas  + f'vvv_xy_f20_c{chip}.txt',xy_vvv_unique,fmt = '%.8f')
        # sys.exit(154)
        # fig, ax = plt.subplots(1,1)
        # ax.scatter(xy_vvv[:,0],xy_vvv[:,1], s =20)
        # ax.scatter(xy_vvv_unique[:,0],xy_vvv_unique[:,1], s =1)
        
        # sys.exit(154)
        p, (pos_img, pos_img_t) = aa.find_transform(xy_s, xy_vvv_unique, max_control_points=200)
        
        print(f'Scale from aa = {p.scale: .2f}')
        print(f'Rotation from aa = {p.rotation*180/np.pi: .2f}º')
        
        
    
        x_overlap = np.array(vvv_overlap['x'][idx])
        y_overlap = np.array(vvv_overlap['y'][idx])
        
    
        indices = []
        
    
        for coord in pos_img_t:
            x, y = coord
            # Find the indices where both x and y match in vvv_overlap
            index = np.where((x_overlap == x) & (y_overlap == y))[0]
            
            # If a match is found, append the index
            if index.size > 0:
                indices.append(index[0])  # Assuming one-to-one match
            else:
                indices.append(-1)  # Append -1 for unmatched points
        
        # Convert indices list to numpy array for easier handling
        indices = np.array(indices)
        
        # Print or return the indices
        print('Comon stars in the C%s = %s'%(chip,len(indices)))
        vvv_com = vvv_overlap[idx][indices]
        vvv_com.write(pruebas + 'vvv_common_c%s.txt'%(chip), format = 'ascii', overwrite = True)
        
        
        gns_com_xy = Table(pos_img, names = ('x','y'))
        vvv_com_xy = Table(pos_img_t, names = ('x','y'))
        gns_com_xy.write(pruebas + 'gns_c%s_com.txt'%(chip), format = 'ascii',overwrite=True)
        vvv_com_xy.write(pruebas + 'vvv_c%s_com.txt'%(chip), format = 'ascii',overwrite=True)
        
    
        vvv_com_ad = SkyCoord(ra = vvv_com['ra'],dec = vvv_com['dec'],unit = 'degree')
        xy_com = np.vstack((pos_img[:,0],pos_img[:,1]))
        wcs_new = fit_wcs_from_points(xy_com, vvv_com_ad, projection="TAN")
        
        data_cube = cube[0].data
        header = cube[0].header  # Get the header from the corresponding extension
        m_cube = np.stack([data_m]*cube[0].header['NAXIS3'], axis = 0)
        for i in range(cube[0].header['NAXIS3']-1):
            image_i +=1
            # Get the data and header for each image (extension)
            if chip ==1:
                data = data_cube[i][0:768,0:2048]# Use i+1 since the first is the PrimaryHDU
                wcs_header = wcs_new.to_header()                
                m_data = m_cube[i][0:768,0:2048]
               
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value
                header['FILTER'] = 'H' # this keyword is needed for Scamp.
                fits.writeto(slices_aligned + '%s_image_%04d.%02d.fits'%(field,image_i,chip),
                             data = data, header = header, overwrite= True)
                
                fits.writeto(slices_aligned + '%s_image_%04d.%02d.weight.fits'%(field,image_i,chip),
                             data = m_data, header = header, overwrite= True)
            if chip == 2:
                data = data_cube[i][0:768,2048:]
                wcs_header = wcs_new.to_header()
                m_data = m_cube[i][0:768,2048:]
                
             
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value
                header['FILTER'] = 'H'        

                fits.writeto(slices_aligned + '%s_image_%04d.%02d.fits'%(field,image_i,chip),
                             data = data, header = header, overwrite= True)
                
                fits.writeto(slices_aligned + '%s_image_%04d.%02d.weight.fits'%(field,image_i,chip),
                             data = m_data, header = header, overwrite= True)
           
            if chip == 3:
                data = data_cube[i][768:,2048:]
                wcs_header = wcs_new.to_header()
                
                m_data = m_cube[i][768:,2048:]
                
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value
                
                header['FILTER'] = 'H'        

                fits.writeto(slices_aligned + '%s_image_%04d.%02d.fits'%(field,image_i,chip),
                             data = data, header = header, overwrite= True)
                
                fits.writeto(slices_aligned + '%s_image_%04d.%02d.weight.fits'%(field,image_i,chip),
                             data = m_data, header = header, overwrite= True)
              
               
            if chip == 4:
                data = data_cube[i][768:,0:2048]
                wcs_header = wcs_new.to_header()
                
                m_data = m_cube[i][768:,0:2048]
                
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value
                
                header['FILTER'] = 'H'        
                
                
                fits.writeto(slices_aligned + '%s_image_%04d.%02d.fits'%(field,image_i,chip),
                             data = data, header = header, overwrite= True)
                
                fits.writeto(slices_aligned + '%s_image_%04d.%02d.weight.fits'%(field,image_i,chip),
                             data = m_data, header = header, overwrite= True)
        last_id = image_i
        image_i -=  cube[0].header['NAXIS3']-1
        
    print(30*'+',f'\nAfter cube{li}, index ={last_id} \n',30*'+')
    dic_sl['l%s'%(li)] = cube[0].header['NAXIS3']-1
    
    
# =============================================================================
#     os.remove(clean + 'cube%s.fits'%(li+1))
# =============================================================================
   
    # if li == 1:
    #     break    
for im in range(1,last_id+1):
        command = f'missfits -SAVE_TYPE REPLACE -OUTFILE_TYPE MULTI -WRITE_XML N -NEXTENSIONS_MIN 4 -SPLIT_SUFFIX .%02d.weight.fits {field}_image_{im:04d}'
        result = subprocess.run(command, check=True, cwd = slices_aligned,shell = True)
        
        rn_c =  f'mv {field}_image_{im:04d}.fits {field}_image_{im:04d}.weight.fits'
        res_rn = subprocess.run(rn_c, check=True, cwd = slices_aligned,shell = True)
                

for im in range(1,last_id+1):
        command = f'missfits -SAVE_TYPE REPLACE -OUTFILE_TYPE MULTI -WRITE_XML N -NEXTENSIONS_MIN 4 -SPLIT_SUFFIX .%02d.fits {field}_image_{im:04d}'
        result = subprocess.run(command, check=True, cwd = slices_aligned,shell = True)


   




        
