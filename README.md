# SeptinManuscriptData
Code for Gabbert, et. al. Septins regulate border cell surface geometry, shape, and motility downstream of Rho in _Drosophila_. 
Developmental Cell 2023.
This code was used to perform 2D geomstats shape analysis, generate 3D surface texture models, and perform spectral decomposition analysis.

# STAR Protocol Code
Some of this code can be used and modified to generate 3D models and perform spectral decomposition analysis with new data, described in detail in the STAR Protocol manuscript.

## Pre-requisites to follow STAR Protocol
- [MATLAB](https://www.mathworks.com/)
  - Deep Learning Toolbox
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Curve Fitting Toolbox   
- [ImageJ/FIJI](https://fiji.sc/) 
- [Ilastik](https://www.ilastik.org/)
- [MeshLab](https://www.meshlab.net/)
- Github scripts and folders
  - STAR_Methods_Tissue_Cartography_Modeling (to generate 3D models)
  - SpectralAnalysis folder
      - script_spectralAnalysis_basic (to perform spectral analysis on a single sample)
      - script_spectralAnalysis_acrossConditions (to compare multiple genotypes or conditions)
  - ImSAnE folder
  - gptoolbox folder
  - external folder

