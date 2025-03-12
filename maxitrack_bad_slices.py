#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:18:34 2024

@author: amartinez
"""

import subprocess
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import os
import re

import csv
part = 1
chip = 1        
prior= 0.05
morralla = '/Users/amartinez/Desktop/morralla/'

all_ext_values = []

for chip in range(1,5):
    folder = '/Volumes/teabag-alvaro/gns_gd/gns1/pruebas/part%s_B1_B6/chip%s/'%(part,chip)
    maxi = Table.read(folder + 'maxitrack_c%s_p%s.out'%(chip, prior), names = ('file','prob'), format = 'ascii')
    
    maxi['file'] = [os.path.basename(path) for path in maxi['file']]
    
    def extract_number(filename):
        match = re.search(r'(\d+)\.miss.fits', filename)
        return match.group(1) if match else None
    
    maxi['ext'] = [extract_number(os.path.basename(path)) for path in maxi['file']]
    
    bad_sl = maxi['prob'] > 0.7
    maxi_bad = maxi[bad_sl]
    maxi_bad['prob'] = [f"{p:.8f}" for p in maxi_bad['prob']]
    
    # Append 'ext' values to the list
    all_ext_values.extend(maxi_bad['ext'])
    


    maxi_bad['prob','ext'].write(morralla + f'bad_sl_c{chip}_p{prior}.csv', overwrite=True)
    np.savetxt(morralla + f'bad_sl_c{chip}_p{prior}_part{part}.csv',np.c_[maxi_bad['prob'], maxi_bad['ext']], fmt='%s',delimiter= ",",header = f'chip{chip} prior{prior} ,ext')
    
    
    # maxi_bad.write(morralla + f'bad_sl_c{chip}_p{prior}.csv', overwrite = True)
# %%
# Convert to integers and remove duplicates
all_ext_values = [int(value) for value in all_ext_values if value is not None]


unique_ext_values, counts =  np.unique(all_ext_values, return_counts = True)

# Print or use the unique 'ext' values
# Save the unique list as a CSV file
# =============================================================================
# csv_file_path = morralla + 'unique_ext_values.csv'
# with open(csv_file_path, mode='w', newline='') as file:
#     writer = csv.writer(file)
#     writer.writerow(['ext'])  # Write the header
#     for value in unique_ext_values:
#         writer.writerow([value])
# =============================================================================

np.savetxt(morralla + f'bad_slices_prior%s_part{part}.csv'%(prior), np.c_[unique_ext_values, counts], delimiter= ",", fmt = "%d", header = f'prior{prior},ext')

# print(f"Unique 'ext' values saved to {csv_file_path}")
# %%

# t = Table()

# for i in range(1,5):
#     bad_s = np.genfromtxt(morralla + f'bad_sl_c{i}_p{prior}.csv', delimiter = ",")
#     if i ==1:
#         t['ext'] = bad_s[:,1]
#         t['chips'] = [1]
#     else:
#         for e in t['ext']:
#             if np.isin(e,t['ext']) == True:
                


from collections import defaultdict

ext_c = defaultdict(list)

for i in range(1,5):
    file = np.loadtxt(morralla + f'bad_sl_c{i}_p{prior}_part{part}.csv',delimiter = ",")
    for e in file[:,1]:
        ext_c[e].append(i)

t = Table()
t[f'ext p{prior}'] = np.array(list(ext_c.keys())).astype(int)
# t['chip'] = list(ext_dict.values())
t['chips'] = [','.join(map(str, values)) for values in ext_c.values()]
t.sort(f'ext p{prior}')
t.write(morralla + f'chips_bad_sl_p{prior}_part{part}.csv', format =  'csv', overwrite = True)



# %%


    
    
    
    
    
    
    
    
    



