# Geometric Distortions in GNS


## There is not yet a define pipeline. Here the name and description of some scripts

## 
- Runs `missfits`on cubes to divide them into slices. Then make a .list file containing the paths and names of the slices and feed with it `maxitrack`
- **% missfits -d > default.missfits** generates de default configurataion file. 
- On the configuration file you can set `OUTFILE_TYPE  SPLIT` and `SAVE_TYPE NEW`. This would slice the MEF and conserve the original file, while saving the slices under a new name.
- Make a .list files with then names and tha absolutes paths of the new genareted *miss* slices: ** ls *.miss.fits | xargs realpath > part[1,1]_c[1,2,3,4]_fits.list **

### maxitract_chips.sh 
- Runs maxitract on fits slices for diffentes chips and rename the output file maxitract.output

### maxitraxt_bad_slices.py

- Combines all the maxitrack output files in a single list. This is combinient for comparing the bad slices selected by maxitrack with the bad slices selected by the operator (me)






# OFFDATED

## Steps to Estimate and Apply Geometric Distortions (GD) on GNS2

---


### 0. Prepare Cubes for GD Corrections
- Run `cubes_for_gd_corrections.py` to group all the cleaned cubes (reduced for sky, flat-fielded, and cleaned for cosmic rays) into four FITS files, one for each chip.
- The generated cubes will have the correct format for use with Astromatic software. **Inspect the generated cubes carefully**.

#### 0.1 Optional: Remove Bad Frames
- Use `delete_bad.py` if additional bad slices are detected after inspection. In principle, the cubes used in `cubes_for_gd_corrections.py` should already be cleaned of bad frames.

### ⚠️ Important Warning:
- The Astromatic software is currently **not working properly** on `teatime`, except for `MISSFITS`.  
- For the subsequent steps, **move the cubes generated in the previous steps to `teabag`** to continue processing.  

## Astromatic.py

- Runs `sextractor`, `scamp` and `SWarp`.

### 1. Run Source Extractor
- Execute `source-extractor cube.fits -c default_c[1,2,3,4].sex`. This will generate a **MEF (Multi-Extension FITS)** file with different pointings, with an extension for each image.
- Ensure this cube is created with the `cubes_for_gd_corrections.py` script.
- Edit `default.param` to select the variables needed by SCAMP. A list of all SExtractor variables can be found in `default.sex`.

### 2. Run SCAMP
- Run the command: `scamp sextractor.cat -c scamp.conf`.
- SExtractor will generate a catalog to feed SCAMP, which calculates the geometric distortion. A `header.head` file will be generated.

### 3. Run SWarp
- Run `SWarp cube.fits`. 
- The `.head` file generated by SCAMP is used to modify the original `cube.fits`. Ensure this file is referenced in the configuration.

### 4. Apply Axis Correction
- Run `axis_correction.py` to equalize the axes in the FITS files generated by SWarp, updating the WCS information accordingly.
- Groups the files in cubes of the same size than those in **cleaned** directory.

---

Follow these steps to properly apply geometric distortion corrections to GNS2 data.


