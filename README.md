# timmw_and_zeddepth
post processing raw mmWave data and depth data

## mmWave SAR data
* **postproc_sar_030623.m**
  * process SAR data to generate 3D intensity maps
  * check if paths to load .bin files are correct
  * if have more than one view, add an outer loop
  
## Depth image data
As of Jan. 2022, depth images are not needed to be cropped before feeding into nueral network. Only need to scale the color to a specific limit (e.g. (255->0.1m, 0->15m)
* **postproc_dep_030623.m**
  * process depth data to generate 2D depth images with color-coded limits
  * the depth data is the one extracted from SVO (using ZED SDK or other latest methods)
