%==========================================================================
% coreg_input()
function R = coreg_input(Nii_x,DoCoReg,ref)
if nargin < 3, ref = [1 1]; end

C = numel(Nii_x);
R = cell(1,C);
for c=1:C
    N       = numel(Nii_x{c}); 
    R{c}    = cell(1,N);
    R{c}(:) = {eye(4)};
end

if numel(size(Nii_x{1}(1).dat)) > 2
    
    Vr = spm_vol(Nii_x{ref(1)}(ref(2)).dat.fname);    
    for c=1:C
        N = numel(Nii_x{c});         
        for n=1:N
            if (c == ref(1) && n == ref(2)) || ~DoCoReg
                continue
            else    
                Vm      = spm_vol(Nii_x{c}(n).dat.fname);
                x       = spm_coreg(Vr,Vm);
                R{c}{n} = spm_matrix(x(:)');
            end        
        end
    end
end
end
%==========================================================================