#!/bin/bash

# Loop over the values 1 and 2
prior=0.05
part=1
for chip in 1 2 3 4
do
    # Define the folder path
    folder="/home/data/alvaro/gns_gd/gns1/pruebas/part${part}_B1_B6/chip${chip}/"

    # Run maxitrack with the specified list file
    maxitrack "${folder}part${part}_c${chip}_fits.list" -v --prior ${prior}
    #maxitrack "${folder}clean_FAST-Bulge-6.Ks.2017-04-05T08_55_22_2.001.miss.fits" -v --prior ${prior}
    
    mv "maxitrack.out" "${folder}maxitrack_c${chip}_p${prior}.out"
    #mv "maxitrack.out" "${folder}maxitrack_c${chip}_DELETE.out"
    

done
 
 
 
