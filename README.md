# A Tool for Super-Resolving Multimodal Clinical MRI

<span style="color:blue">**NEW: Python version (using PyTorch) now available from: https://github.com/brudfors/UniRes**</span>

Multi-channel total variation (MTV) super-resolution of magnetic resonance imaging (MRI) data. The tool can be run either as a Docker image or MATLAB (using SPM12).

The code is based on the algorithm described in the papers:

     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     A Tool for Super-Resolving Multimodal Clinical MRI.
     2019 arXiv preprint arXiv:1909.01140.     
     
     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
     MRI Super-Resolution Using Multi-channel Total Variation.
     In Annual Conference on Medical Image Understanding and Analysis
     2018 Jul 9 (pp. 217-228). Springer, Cham.             

## Using Docker

With Docker ([https://docs.docker.com/get-started/](https://docs.docker.com/get-started/)) installed on your system, use the following steps to super-resolve a set of subject MRI scans to 1 mm isotropic voxel size:

1. Get the spm_superres Docker image: 
``` bash
docker pull mbrud/spm_superres
```

2. Open an spm_superres container: 
``` bash
docker run -ti --rm -v [PTH-IN]:/input mbrud/spm_superres
```  
where `[PTH-IN]` is the full path to a directory containing a set of subject MR images (OBS: all NIfTIs will be read from this folder, so make sure that it only contains images you want super-resolved).

3. Execute the following command in the container: 
``` bash
/opt/spm12/spm12 function spm_superres input
```

4. When the algorithm has finished, you will find the super-resolved images in the `[PTH-IN]` folder, prefixed `'y'`.

**OBS!** This algorithm can end up using A LOT of RAM. If you are running the Docker version on a Mac, the default memory allocation is just 2GB, which can lead to a cryptic *killed* error message in the terminal. Increasing the available docker memory can help:

https://stackoverflow.com/questions/44417159/docker-process-killed-with-cryptic-killed-message

## Using MATLAB

If you want to change model parameters (e.g., the super-resolved images' voxel size, or increase the regularisation) or modify the code itself you will need to run it using MATLAB. If so, then just pull/download the code from this repository onto your computer. Make sure that the SPM12 software is on Matlab's path. The most recent version of SPM12 can be downloaded from [www.fil.ion.ucl.ac.uk/spm](http://www.fil.ion.ucl.ac.uk/spm/). If you get error messages when running the code, it is probably because your SPM version is too old. 

The following is an example of super-resolving three MR scans of a subject to 1 mm isotropic voxel size:
``` matlab
% Paths to some MR images of the same patient (in nifti format)
P    = cell(1,3);
P{1} = 'MRI-contrast1.nii';
P{2} = 'MRI-contrast2.nii';
P{3} = 'MRI-contrast3.nii';

spm_superres(P);
```
Output images are written to the same folder as the input images, prefixed `'y'`.

### Improved runtime (Linux and Mac)

For a faster algorithm, consider compiling SPM with OpenMP support. Just go to the *src* folder of SPM and do:
``` bash
make distclean
make USE_OPENMP=1 && make install
```

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.
