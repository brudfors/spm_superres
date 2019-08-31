%==========================================================================
% preproc_input()
function Nii = preproc_input(Nii,DirOut,NumWorkers,ref,prefix)
if nargin < 4, ref    = [1 1]; end
if nargin < 5, prefix = 'r';   end

if ~(exist(DirOut,'dir') == 7) 
    mkdir(DirOut);
end
    
C = numel(Nii);

%---------------------------
% Copy to DirOut
%---------------------------

for c=1:C
    N = numel(Nii{c});
    for n=1:N            
        fname       = Nii{c}(n).dat.fname;
        [~,nam,ext] = fileparts(fname);
        nfname      = fullfile(DirOut,[nam ext]);
        copyfile(fname,nfname);
        Nii{c}(n)   = nifti(nfname);
    end
end

%---------------------------
% Coregister and reslice to same size
%---------------------------

reffname = Nii{ref(1)}(ref(2)).dat.fname;

parfor (c=1:C,NumWorkers)
    N = numel(Nii{c});
    for n=1:N            
        if c == ref(1) && n == ref(2), continue; end
            
        fname         = Nii{c}(n).dat.fname;
        [pth,nam,ext] = fileparts(fname);

        matlabbatch                                                 = {};
        matlabbatch{1}.spm.spatial.coreg.estwrite.ref               = {reffname};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source            = {fname};      
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep      = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm     = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp   = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap     = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask     = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix   = prefix;          
        
        spm_jobman('run',matlabbatch);
        
        nfname    = fullfile(pth,[prefix nam ext]);
        Nii{c}(n) = nifti(nfname);
        delete(fname);
    end
end
end
%==========================================================================