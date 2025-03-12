from astropy.io import fits
import os
import sys
import glob



# Path where the FITS files are stored


for chip in range(3,5):
    path = '/home/data/alvaro/gns_gd/gns1/H/B6/cubes_gd/slices/chip%s/'%(chip)
    n_file = glob.glob(path + '*.weight.fits')
    print(len(n_file))
    
    # Loop over the series of FITS files
    for i in range(1, len(n_file)+1):  # adjust the range based on how many files you have
        filename = f"B6_image_c{chip}.{i:04d}.weight.fits"  # Assuming 4 digits (0001, 0002, ..., 0020)
        filepath = os.path.join(path, filename)
        
        try:
            # Open the FITS file
            with fits.open(filepath, mode='update') as hdul:
                # Add the 'FILTER' keyword with value 'H' to the header of the primary HDU
                hdr = hdul[0].header
                hdr['FILTER'] = 'H'
                print(f"Added FILTER keyword to {filename}")
                
                # Save changes
                hdul.flush()
        
        except Exception as e:
            print(f"Error processing {filename}: {e}")
        
    
    # Loop over the series of FITS files
    for i in range(1, len(n_file)+1):  # adjust the range based on how many files you have
        filename = f"B6_image_c{chip}.{i:04d}.fits"  # Assuming 4 digits (0001, 0002, ..., 0020)
        filepath = os.path.join(path, filename)
        
        try:
            # Open the FITS file
            with fits.open(filepath, mode='update') as hdul:
                # Add the 'FILTER' keyword with value 'H' to the header of the primary HDU
                hdr = hdul[0].header
                hdr['FILTER'] = 'H'
                print(f"Added FILTER keyword to {filename}")
                
                # Save changes
                hdul.flush()
        
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            
