# A Tool for Super-Resolving Multimodal Clinical MRI

Multi-channel total variation (MTV) super-resolution of magnetic resoinace imaging (MRI) data. The tool can be run either as a Docker image or MATLAB (using SPM12).

The code is based on the algorithm described in the paper:

     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     MRI Super-Resolution Using Multi-channel Total Variation.
     In Annual Conference on Medical Image Understanding and Analysis
     2018 Jul 9 (pp. 217-228). Springer, Cham.   

## Using Docker

Coming soon.

## Using MATLAB

Just pull/download the code from this repository onto your computer. Make sure that the SPM12 software is on Matlab's path. The most recent version of SPM12 can be downloaded from [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/). If you get error messages when running the code, it is probably because your SPM version is too old. 

The following is an example of super-resolving three MR scans of a subject to 1 mm isotropic voxel size:
~~~~
% Paths to some MR images of the same patient (in nifti format)
P    = cell(1,3);
P{1} = 'MRI-contrast1.nii';
P{2} = 'MRI-contrast2.nii';
P{3} = 'MRI-contrast3.nii';

spm_superres(P);
~~~~
Output images are written to the same folder as the input images, prefixed `'y'`.

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.
