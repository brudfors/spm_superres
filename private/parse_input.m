%==========================================================================
% parse_input()
function Nii = parse_input(pths)
C   = numel(pths);
Nii = cell(1,C);
for c=1:C
    Nii{c} = nifti(pths{c});
end
end
%==========================================================================