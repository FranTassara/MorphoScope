Synthetic Volumes for MorphoScope Validation
=============================================

Generated volumes to demonstrate and validate the spread quantification algorithm.

Volume Specifications
---------------------
- XY dimensions: 800 x 800 pixels
- Z dimensions: 20 slices
- Pixel size (XY): 0.1 um
- Z step: 1.0 um
- Physical size: 80.0 x 80.0 x 20.0 um

Demo Categories
---------------

1. spread_concept/
   Demonstrates what spread measures:
   - Spheres of different sizes: spread scales with size
   - Ellipsoids elongated in different axes: anisotropic spread
   - Same geometry at different intensities: spread unchanged
   
2. pca_rotation/
   Demonstrates PCA axis alignment:
   - Same ellipsoid at different rotations
   - After PCA alignment, all should give similar spread values
   
3. fasciculation/
   Demonstrates the biological phenomenon:
   - fasciculated (dusk state): axons bundled tightly, low Y/Z spread
   - defasciculated (dawn state): axons spread out, high Y/Z spread
   
4. axonal_volume/
   Demonstrates axonal volume measurement:
   - Same size, different intensity: different volume
   - Different sizes: different volume
   - Shell vs solid: spread similar, volume different

Usage with MorphoScope
----------------------
1. Load any .tif file from these directories
2. The entire image can be used as ROI (no selection needed)
3. Compare spread values between related volumes
4. Use to verify algorithm correctness against known geometries

Expected Results
----------------
spread_concept/:
  - Larger spheres: larger spread values
  - X-elongated ellipsoid: largest X spread
  - Different intensity: same spread (intensity-weighted mean)

pca_rotation/:
  - All rotations should give similar spread after PCA alignment
  - X spread should always be largest (PCA aligns to maximum variance)

fasciculation/:
  - Fasciculated: low Y/Z spread, similar X spread
  - Defasciculated: high Y/Z spread
  - Branching structures: intermediate values

axonal_volume/:
  - Higher intensity: higher axonal volume
  - Larger structure: higher axonal volume
  - Shell vs solid: similar spread, different volume


Real Volumes for MorphoScope Validation
=============================================

MorphoScope has been validated against the original MATLAB implementation 
(Petsakou et al., 2015):

| Metric | Mean difference | Pearson r |
|--------|-----------------|-----------|
| Spread X | 0.4% | 0.9997 |
| Spread Y | 3.0% | 0.9911 |
| Spread Z | 1.4% | 0.9982 |

Validation performed using real confocal data (n=4 samples, ZT2 and ZT14).
