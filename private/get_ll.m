%==========================================================================
% get_ll()
function [ll1,ll2,ll] = get_ll(Nii_x,Nii_y,Nii_Dy,tau,dat,NumWorkers,Inplane1mm)
% log posterior
C   = numel(Nii_y);
ll1 = 0;
ll2 = 0;
parfor (c=1:C,NumWorkers)
% for c=1:C, fprintf('OBS: for!\n')
    y = get_nii(Nii_y{c});
    for n=1:numel(Nii_x{c})
        datcn = dat(c).A(n);
        x     = get_nii(Nii_x{c}(n));                                
        
        if Inplane1mm && datcn.do_pm
            % Downsample observed data to have 1 mm inplane
            % resolution
            y0 = double(apply_affine(datcn.D,datcn.dm(1:3)));
            x  = spm_bsplins(double(x),y0(:,:,:,1),y0(:,:,:,2),y0(:,:,:,3),[0 0 0 0 0 0]);   
            x  = single(x);
            y0 = [];
        end      
                
        msk = get_msk(x);
        Ay  = pm('A',y,dat(c).A(n));
        ll1 = ll1 + (tau{c}(n)/2)*sum(sum(sum((double(x(msk)) - double(Ay(msk))).^2)));
        x   = [];
        Ay  = [];
    end
    y = [];    
    
    Dy  = get_nii(Nii_Dy{c});
    ll2 = ll2 + sum(Dy.^2,4);
    Dy  = [];
end
ll2 = sum(sum(sum(double(sqrt(ll2)),3),2),1);
ll1 = -ll1;
ll2 = -ll2;
ll  = ll1 + ll2;
end
%==========================================================================