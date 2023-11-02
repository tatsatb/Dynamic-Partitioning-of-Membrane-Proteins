# Dynamic Partitioning of Membrane Proteins

This repository hosts the code that was developed for the article: *“A dynamic partitioning mechanism polarizes membrane protein distribution”* (DOI: [10.1101/2023.01.03.522496](https://doi.org/10.1101/2023.01.03.522496)). Please see the *Methods* section of the paper for details on development of the computational modeling and image analysis workflow. 

This repository contains two main sections:

- The folder `Computational_Model` contains code to perform two-dimensional spatiotemporal simulation using an excitable network-based stochastic reaction-diffusion system. These programs use [URDME](https://github.com/URDME/urdme) framework. The code is this folder was used to generate (the modified verison of) Figure 6, Supplementary Figure 9, and Supplementary figure 10 of the above mentioned paper. 

- The folder `Image_Analysis` contains code for computing intensity profiles and performing optical flow analysis that was developed to quantitate photoconversion microscopy-based protein tracking assay. This code was used to generate Figure 4, Supplementary Figure 6, and Supplementary figure 7 of the above mentioned paper. 

Please read the [citation](#Citation/Restrictions) details below if you want to use/incorporate/modify a part of this repository in your research. 

## File details

Please go through the readme files in `Computational_Model` ([here](/Computational_Model/ComputationalModel.md)) and `Image_Analysis` ([here](/Image_Analysis/ImageAnalysis.md)) folders to learn about system requirements, dependencies, directory details, and quickstart instructions. 


## Authors

The code is this repository was developed by Debojyoti Biswas, Pablo A. Iglesias, and Tatsat Banerjee (Johns Hopkins University, Baltimore, MD, USA). 

## Citation/Restrictions

This program is a free software (please see [License](#license) for details). However, if you are using this code in your work, please cite our work as:


> **Tatsat Banerjee, Satomi Matsuoka, Debojyoti Biswas, Yuchuan Miao, Dhiman Sankar Pal, Yoichiro Kamimura, Masahiro Ueda, Peter N. Devreotes, and Pablo A. Iglesias** _“A dynamic partitioning mechanism polarizes membrane protein distribution”_, bioRxiv, 2023. DOI: [10.1101/2023.01.03.522496](https://doi.org/10.1101/2023.01.03.522496).

_Note:_ Please visit the [GitHub repository ](https://github.com/tatsatb/Dynamic-Partitioning-of-Membrane-Proteins) for updated citation details. 

## License 

Copyright © 2023 Tatsat Banerjee, Satomi Matsuoka, Debojyoti Biswas, Yuchuan Miao, Dhiman Sankar Pal, Yoichiro Kamimura, Masahiro Ueda, Peter N. Devreotes, and Pablo A. Iglesias.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 