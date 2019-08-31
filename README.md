# spm_superres

Multi-channel total variation (MTV) super-resolution of MR data using SPM12.

The code is based on the method described in the paper:

     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     MRI Super-Resolution Using Multi-channel Total Variation.
     In Annual Conference on Medical Image Understanding and Analysis
     2018 Jul 9 (pp. 217-228). Springer, Cham.
     
The most up-to-date PDF version of the paper is available from https://arxiv.org/abs/1810.03422.

## Dependencies

This project has strong dependencies on SPM12 and its `Shoot` and `Longitudinal` toolboxes. These should all be added to Matlab's path. The most recent version of SPM can be downloaded from [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/). If you get error messages when running the code, it is probably because your SPM version is too old.

## Example

~~~~
addpath('spm_superres')

% Paths to some MR images of the same patient (in nifti format)
P    = cell(1,3);
P{1} = 'MRimage1.nii';
P{2} = 'MRimage2.nii';
P{3} = 'MRimage3.nii';

spm_superres(P);
~~~~
