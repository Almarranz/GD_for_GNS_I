# Geometric Distortions in GNS

## Steps to Estimate and Apply Geometric Distortions (GD) on GNS2

---

### 0. Prepare Cubes for GD Corrections
- Run `cubes_for_gd_corrections.py` to group all the cleaned cubes (reduced for sky, flat-fielded, and cleaned for cosmic rays) into four FITS files, one for each chip.
- The generated cubes will have the correct format for use with Astromatic software. **Inspect the generated cubes carefully**.

#### 0.1 Optional: Remove Bad Frames
- Use `delete_bad.py` if additional bad slices are detected after inspection. In principle, the cubes used in `cubes_for_gd_corrections.py` should already be cleaned of bad frames.

### WARNING
____
- For the moment, the Astromatic sofware is not woriking properly in `teatime`, with the exception of `MISSFITS`. 
- For this point on you have to move the cubes geneated in the previous steps to `**teabag**`
___



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

---

Follow these steps to properly apply geometric distortion corrections to GNS2 data.
