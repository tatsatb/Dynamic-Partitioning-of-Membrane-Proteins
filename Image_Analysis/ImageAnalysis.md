## Requirements

The MATLAB 2021a (or later versions) should be installed in the system to execute the programs.  

Additionally, the following MATLAB toolboxes are required: 

1. *Image Processing Toolbox*	(version 11.3 or later)
2. *Computer Vision Toolbox*	(version 10.0 or later).

Please note that these programs were developed in a Windows 10 environment and hence these may not execute correctly in Mac OS or Linux systems. 

## File and Directory Details 

- `Image_Sample` is the primary folder that contains all the required programs and a sample image (`Image_Sample.tif`).

- `Image_Sample.tif` file is provided as an example. Any other images can be used as replacement (it should be put in the root folder).

- `Optical_flow.m` is the main program. It calls other auxiliary functions to conduct the required preprocessing, segmentation, intensity calculations, and optical flow analysis. 

- `Save_Folder` should exist in the root directory so that mask figures, quiver figures, areas of interests, and final data structures can be saved under proper subdirectory. 
