%==========================================================================
% get_nii()
function img = get_nii(nii)
if isa(nii, 'nifti')
    img = single(nii.dat());
elseif isstruct(nii)
    if isfield(nii, 'private')
        img = single(spm_read_vols(nii));
    else
        img = single(nii.dat());
    end
else
    img = single(nii());
end
%==========================================================================