%==========================================================================
% put_nii()
function nii = put_nii(nii,img)
if isa(nii, 'nifti')
    nii.dat(:) = img(:);
elseif isstruct(nii)
    if isfield(nii, 'private')
        nii = spm_write_vol(nii,img);
    else
        nii.dat(:) = img(:);
    end
else
    nii = img;
end
%==========================================================================