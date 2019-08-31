%==========================================================================
% max_bb_orient()
function [mat,dm,vx] = max_bb_orient(Nii,R,D,vx)

if nargin < 4, vx = [1 1 1]; end

if numel(vx) == 1, vx = vx*ones(1,3); end
    
mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for c=1:numel(Nii)
    N = numel(Nii{c});
    for n=1:N
        dmn = size(Nii{c}(n).dat);
        
        if numel(dmn) == 2
            dmn(3) = 0; 
        end

        dmn  = floor(D{c}{n}(1:3,1:3)*dmn(1:3)')';

        t = uint8(0:7);
        s = diag(dmn+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                              bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                              bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
                          
        mat = R{c}{n}\(Nii{c}(n).mat/D{c}{n});
        
        s  = bsxfun(@plus,mat(1:3,1:3)*s,mat(1:3,4));
        mx = max(mx,max(s,[],2));
        mn = min(mn,min(s,[],2));
    end
end

mat = spm_matrix(mn-1)*diag([vx 1])*spm_matrix(-[1 1 1]);
dm  = ceil((mat\[mx'+1 1]')');
dm  = dm(1:3);

if dmn(3) == 0
    dm(3) = 1;
end
end
%==========================================================================