## Requirements

The codes in this folder was developed with MATLAB 2022a on macOS (version 12) which contains Xcode. 

The programs have been tested with: 

MATLAB 2021a on a Linux Mint 20.2 Uma Cinnamon Edition (based on Ubuntu 20.04 LTS) environment which had pre-built GNU Compiler Collection (GCC) 9.4.0. 

Other Linux distros and versions may be compatible; however, please note that a supported C compiler [[Details on MathWorks' webpage](https://www.mathworks.com/support/requirements/previous-releases.html)] need to be installed on the system that MATLAB can access to compile C source codes to generate MEX files. The codes in this folder uses [URDME](https://github.com/URDME/urdme) framework and hence, Microsoft Windows operating systems are _not_ supported. 

Additionally, the following toolboxes/products of MATLAB are required dependencies:  

1. Image Processing Toolbox	(version 11.3 or later)
2. Partial Differential Equation Toolbox (version 3.6 or later)
3. Parallel Computing Toolbox' (version	7.4 or later)
4. MATLAB Parallel Server (version 7.4 or later)
5. Polyspace Bug Finder (version 3.4 or later).

## File and Directory Details

- `Partitioning_Simulation.m` is the file that perform the simulation of an excitable network by incorporating experimentally obtained differential diffusion profiles of lipid-anchored and peripheral membrane proteins inside different states of the membrane.   

- `Control_Simulation.m` is the file that perform simulation of an excitable network but uses same diffusion coefficients for bound and unbound forms of membrane proteins (instead of using experimentally obtained values).

- `Initial_Condition_Data.mat` loads the initial conditions for the simulation (for either simulation). 

- `urdme_local` folder contains URDME package. It needs to be in the PATH for running `Partitioning_Simulation.m` or `Control_Simulation.m` file.  