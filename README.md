# enhanceIncrementalTeethLines
MATLAB software to automatically enhances incremental teeth lines on X-ray microtomography images. The method is described in Tanner et al., "Extended-field synchrotron microtomography for non-destructive analysis of incremental lines in archeological human teeth cementum", Proceedings of SPIE 11840 (2021) 1184019, DOI: [10.1117/12.2595180](https://doi.org/10.1117/12.2595180).

The scripts are written in Matlab.

# Dataset
Please find an exemplary microtomography dataset (T1, height step 2, tooth4_hs02) in the zenodo repository DOI: [10.5281/zenodo.11029099](https://doi.org/10.5281/zenodo.11029099). It should be unzipped and the tif files reco_????.tif placed in dataDir1/tooth4_hs02/reco/ before running the programs.

# Usage
Define data directories and parameters in dataParameterDefinitionTeeth.m, e.g. set sampleIdx=1 and hs=2. Run extractCementumPatchesEnhanceIL.m to process the data.

# Main programs
data directories and parameters are defined in dataParameterDefinitionTeeth.m.

extractCementumPatchesEnhanceIL.m extracts tooth cementum patches and enhances incremental line appearance automatically using the follow main steps
1) extract patches from the outer band containing the incremental teeth lines from the image of the whole tooth 
2) straighten patches
3) automatic enhancement of incremental lines

optimizeProjection.m does the automatic enhancement of the incremental teeth lines by finding a rotation angle which maximized standard devation of 1D intensity profiles orthogonal to the teeth lines after mean intensity projecton.
