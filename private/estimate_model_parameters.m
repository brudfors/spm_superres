%==========================================================================
% estimate_model_parameters()
function [tau,lam,rho] = estimate_model_parameters(Nii,LamScl,RhoScl,NumWorkers,Verbose,LenaNoiseStd)
if nargin < 6, LenaNoiseStd = 1e1; end   % For 'lena', std of additive Gaussian noise

C   = numel(Nii);
tau = cell(1,C); 
lam = cell(1,C);
parfor (c=1:C,NumWorkers)
    N        = numel(Nii{c});
    NoiseStd = [];
    mu       = [];
    if ~isa(Nii{c}(1),'nifti')
        % Lena image
        n        = 1;
        x        = get_nii(Nii{c}(n));
        NoiseStd = LenaNoiseStd;
        tau{c}   = 1/(NoiseStd.^2);            
        mu       = double(mean(x(isfinite(x)))); 
    else
        % MRI stored in nifti
        for n=1:N
            [NoiseStd,mu] = spm_noise_estimate(Nii{c});
            tau{c}        = 1./(NoiseStd.^2);          
        end
    end
    lam{c} = LamScl/double(mean(mu));
    
    % Scale with number of observations to give a little more weight to the
    % prior when there are multiple observations of the same channel
    lam{c} = sqrt(N)*lam{c}; 
    
    if Verbose
        for n=1:N
            fprintf('c=%i, i=%i | sd=%f, mu=%f -> tau=%f, lam=%f\n', c, n, NoiseStd(n), mu(n), tau{c}(n), lam{c});
        end
    end
end

atau = [];
for c=1:C
    N = numel(tau{c});
    for n=1:N
        atau = [atau tau{c}(n)];
    end
end

% This value of rho seems to lead to reasonably good convergence
rho = RhoScl*sqrt(mean(atau(:)))/mean(cell2mat(lam));
if Verbose
    fprintf('step-size -> rho=%f\n', rho);
end
end
%==========================================================================