%==========================================================================
% get_dat()
function [dat,dm,mat,vx] = get_dat(Nii_x,R,VoxSize,Inplane1mm,Denoise)

C = numel(Nii_x);
D = cell(1,C);
for c=1:C
    N       = numel(Nii_x{c}); 
    D{c}    = cell(1,N);
    D{c}(:) = {eye(4)};
end
if Inplane1mm && ~Denoise
    %---------------------------
    % Make sure that inplane voxel size is one
    %---------------------------
    
    D = cell(1,C);
    for c=1:C
        N    = numel(Nii_x{c});
        D{c} = cell(1,N);
        for n=1:N
            mat = Nii_x{c}(n).mat;
            vxc = sqrt(sum(mat(1:3,1:3).^2));

            % Get down-sampling factor
            d = vxc(1:2);
            if d(1)>=1, d(1) = 1; end
            if d(2)>=1, d(2) = 1; end
            d(3) = 1;

            % NN downsampling    
            D1      = diag([d, 1]);   
            D{c}{n} = D1;
        end
    end                             
end
    
%---------------------------
% Check if projection matrices need to be used
%---------------------------

vx = [];
dm = [];
for c=1:C
    N = numel(Nii_x{c});
    for n=1:N
        vxc = sqrt(sum(Nii_x{c}(n).mat(1:3,1:3).^2)); % Voxel sizes   
        dmc = Nii_x{c}(n).dim;                        % Image size
        dmc = [dmc 1];

        vx = [vx; round(vxc,2)];
        dm = [dm; dmc];
    end
end

if VoxSize == 0
    % Set super-resolved images voxel size to smallest available
    VoxSize = min(min(vx));
end

vx = unique(vx,'rows');
dm = unique(dm,'rows');

%---------------------------
% Build projection matrix struct (dat)
%---------------------------

use_projmat = ~(size(dm,1) == 1 && all(vx == 1)) && ~Denoise;

if ~use_projmat
    % No need to use projection matrices    
    dat = struct;
    for c=1:C
        N = numel(Nii_x{c});
        for n=1:N
            dat(c).A(n).do_pm = false;
        end
    end
    mat = Nii_x{1}(1).mat; % Orientation matrix are the same for all images
else
    % Images are different size and/or different voxel size, use projection
    % matrices         
    [mat,dm,vx] = max_bb_orient(Nii_x,R,D,VoxSize);    
    dat         = init_dat(Nii_x,mat,dm,R,D);    
end
end
%==========================================================================