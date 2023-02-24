function oNii = spm_superres(pths,opt)
% Multi-channel total variation (MTV) super-resolution of MR data.
%
% Requires that the SPM software is on the MATLAB path.
% SPM is available from: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
%
% FORMAT oNii = spm_superres(pths,opt)
%
% INPUT
% pths  - [C x 1] cell array, where each element of the cell array holds a
%                 char array of I_c paths to nifti files (observations of 
%                 channel c). For example, if you for one patient have 
%                 2xT1w, 3xT2w and 1xPDw images, then pths should be 
%                 defined as:
%                         pths{1} = char('T1w1.nii','T1w2.nii')
%                         pths{2} = char('T2w1.nii','T2w2.nii','T2w3.nii')
%                         pths{3} = char('PDw1.nii')
%        or
%                 The path to a folder that can have two structures:
%                 1. The folder contains C niftis of different channels.
%                 1. Each subfolder of the folder contains a set of niftis 
%                    of the same channel. So there are as many subfolders
%                    as channels.
% opt   - Algorithm options, more info below.
%
% OUTPUT
% oNii - nifti object of C processed MRIs
%__________________________________________________________________________
% OPTIONS
% LamScl      - Scaling of regularisation parameter (lambda)           [10]
% RhoScl      - Scaling of step-size parameter (rho)                    [1]
% MaxNiter    - Max number of iterations                               [30]
% NiterNewton - Newton iterations                                       [1]
% Tolerance   - Convergence threshold                                [1e-4]
% DoMTV       - Run either MTV or indepentend TV denoising           [true]
% DirOut      - Directory where to write output (and temporary files, 
%               which are deleted at end of algorithm)                 ['']
% Nii_y0      - Clean reference image                                    []
% NumWorkers  - Number of parfor workers                                [8]
% Verbose     - Show stuff                                              [1]
% DoCoReg     - Do preprocessing (register)                          [true]
% ShowZoomed  - Show one image, zoomed in                           [false]
% MaxMem      - Memory limit to allocate variables as niftis         [4096]
% VoxSize     - Reconstruction voxel size                               [1]
%               If 0, set to smallest available
% Inplane1mm  - Downsample inplane resolution to 1 mm                [true]
% Denoise     - Do just denoising, without super-resolving          [false]
%__________________________________________________________________________
% The general principles are described in the following paper:
%
%     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
%     A Tool for Super-Resolving Multimodal Clinical MRI.
%     2019 arXiv preprint arXiv:1909.01140.   
%
%     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
%     MRI Super-Resolution Using Multi-channel Total Variation.
%     In Annual Conference on Medical Image Understanding and Analysis
%     2018 Jul 9 (pp. 217-228). Springer, Cham.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 2, opt  = struct; end

% Add required SPM toolboxes to path
pth = fileparts(which('spm'));
if ~isdeployed, addpath(pth); end
if ~isdeployed, addpath(fullfile(pth,'toolbox','Longitudinal')); end
if ~isdeployed, addpath(fullfile(pth,'toolbox','Shoot')); end

% Set up boundary conditions
spm_diffeo('bound',1); % OK?
spm_field('bound',1);  % match the gradient operator

%---------------------------
% Options
%---------------------------

opt = spm_superres_lib('get_opt',opt);

LamScl      = opt.LamScl;
RhoScl      = opt.RhoScl;
MaxNiter    = opt.MaxNiter;
NiterNewton = opt.NiterNewton;
Tolerance   = opt.Tolerance;
DoMTV       = opt.DoMTV;
DirOut      = opt.DirOut;
Nii_y0      = opt.Nii_y0;
NumWorkers  = opt.NumWorkers;
Verbose     = opt.Verbose;
DoCoReg     = opt.DoCoReg;
ShowZoomed  = opt.ShowZoomed;
MaxMem      = opt.MaxMem;
VoxSize     = opt.VoxSize;
Inplane1mm  = opt.Inplane1mm;
Denoise     = opt.Denoise;

% Repeatable random numbers
% rng('default');
% rand(1);

%---------------------------
% Parse input data
%---------------------------
        
Nii_x = spm_superres_lib('parse_input',pths);

%---------------------------
% Estimate rigid alignment matrices
%---------------------------

R = spm_superres_lib('coreg_input',Nii_x,DoCoReg);
    
%---------------------------
% Get projection matrix struct and image properties
%---------------------------

[dat,dm,mat,vx] = spm_superres_lib('get_dat',Nii_x,R,VoxSize,Inplane1mm,Denoise);
    
C  = numel(Nii_x); % Number of channels
N0 = 0;            % Total number of observations
for c=1:C
    N0 = N0 + numel(Nii_x{c});
end
if C == 1
    NumWorkers = 0;
end

% Show input
spm_superres_lib('show_stuff',Nii_y0,'y',1,Verbose,ShowZoomed);
if dat(1).A(1).do_pm
    spm_superres_lib('show_stuff',Nii_x, ['x (C=' num2str(C) ' N=' num2str(N0) ')'],2,Verbose,true);
else
    spm_superres_lib('show_stuff',Nii_x, ['x (C=' num2str(C) ' N=' num2str(N0) ')'],2,Verbose,ShowZoomed);
end

%---------------------------
% Decide if to allocate temp vars in niftis or not
%---------------------------

WriteTmpNii = spm_superres_lib('check_do_write',C,dm,MaxMem,NumWorkers);

%---------------------------
% Estimate model parameters
%---------------------------

[tau,lam,rho] = spm_superres_lib('estimate_model_parameters',Nii_x,LamScl,RhoScl,NumWorkers,Verbose);

%---------------------------
% Init variables (creates a load of niftis that will be deleted at end of
% algorithm)
%---------------------------
    
[Nii_y,Nii_z,Nii_w,Nii_Dy,Nii_H] = spm_superres_lib('alloc_vars',WriteTmpNii,Nii_x,dm,mat,DirOut,Verbose);

%---------------------------
% Run algorithm 
%---------------------------

if Verbose
    if Tolerance == 0
        fprintf('\nRunning %d iterations:\n', MaxNiter);
    else
        fprintf('\nRunning max %d iterations:\n', MaxNiter);
    end
end

tic
ll = -Inf;
for it=1:MaxNiter
   
    %---------------------------
    % Update y
    %---------------------------
        
    parfor (c=1:C,NumWorkers)
%     for c=1:C, fprintf('OBS: for!\n')
        N = numel(Nii_x{c});
        y = [];
        for it1=1:NiterNewton
                    
            y = spm_superres_lib('get_nii',Nii_y{c});
            z = spm_superres_lib('get_nii',Nii_z{c});
            w = spm_superres_lib('get_nii',Nii_w{c});
           
            % Compute gradient and Hessian
            gr = single(0);
            H  = single(0);
            
            % Conditional part (mask)
            for n=1:N                
                
                datcn = dat(c).A(n);
                
                % Mask missing data
                x = spm_superres_lib('get_nii',Nii_x{c}(n));
                if Inplane1mm && datcn.do_pm
                    % Downsample observed data to have 1 mm inplane
                    % resolution
                    y0 = double(spm_superres_lib('apply_affine',datcn.D,datcn.dm(1:3)));
                    x  = spm_bsplins(double(x),y0(:,:,:,1),y0(:,:,:,2),y0(:,:,:,3),[0 0 0 0 0 0]);   
                    x  = single(x);
                    y0 = [];
                end                                
                
                % Gradient part
                msk       = spm_superres_lib('get_msk',x);                
                gr1       = tau{c}(n)*spm_superres_lib('pm','At',(spm_superres_lib('pm','A',y,datcn) - x),datcn);
                gr1(~msk) = 0;
                gr        = gr + gr1;
                x         = [];
                
                % Hessian part
                if it == 1 && NiterNewton == 1
                    % Compute only once (does not change)
                    H1       = tau{c}(n)*spm_superres_lib('pm','At',spm_superres_lib('pm','A',ones(dm,'single'),datcn),datcn);
                    H1(~msk) = 0;               
                    H        = H + H1;
                end
                
            end                        
            gr1   = []; 
            H1    = []; 
            msk   = [];
            datcn = [];
            
            if it == 1 && NiterNewton == 1
                % Save Hessian
                Nii_H{c} = spm_superres_lib('put_nii',Nii_H{c},H);
            else
                % Load Hessian
                H = spm_superres_lib('get_nii',Nii_H{c});
            end
            
            gr = gr + spm_superres_lib('diffoperator',w - rho*z,dm,vx,lam{c},'Dt');
            z  = [];
            w  = [];
            
            gr = gr + spm_field('vel2mom',y,[vx 0 rho*lam{c}^2 0]);

            % Do update with spm_field
            y = y - spm_field(H,gr,[vx 0 rho*lam{c}^2 0  1 1]);
            g = []; 
            H = [];
            
            % Ensure non-negative
            y(y < 0) = 0;
            
            % Update nii
            Nii_y{c} = spm_superres_lib('put_nii',Nii_y{c},y);
        end

        % Compute Dy
        Dy        = spm_superres_lib('diffoperator',y,dm,vx,lam{c},'D');
        Nii_Dy{c} = spm_superres_lib('put_nii',Nii_Dy{c},Dy);
        y         = [];
        Dy        = [];
    end    
    
    spm_superres_lib('show_stuff',Nii_y,'yhat',4,Verbose,ShowZoomed);
    
    %---------------------------
    % Update z
    %---------------------------
    
    if DoMTV
        % MTV
        znorm = single(0);
        parfor (c=1:C,NumWorkers)     
            w     = spm_superres_lib('get_nii',Nii_w{c});
            Dy    = spm_superres_lib('get_nii',Nii_Dy{c});
            znorm = znorm + sum((Dy + w/rho).^2,4);        
        end
        Dy = [];
        w  = [];
        
        znorm = sqrt(znorm);
        mtv   = max(znorm - 1/rho,0)./(znorm + eps);
        znorm = [];

        parfor (c=1:C,NumWorkers)    
            w  = spm_superres_lib('get_nii',Nii_w{c});
            Dy = spm_superres_lib('get_nii',Nii_Dy{c});
            z = bsxfun(@times, mtv, (Dy + w/rho));
            
            % Update nii
            Nii_z{c} = spm_superres_lib('put_nii',Nii_z{c},z);
        end    
        Dy = [];
        w  = [];
        z  = [];
        
        spm_superres_lib('show_stuff',mtv,'mtv',3,Verbose,ShowZoomed);
    else
        % Regular TV
        parfor (c=1:C,NumWorkers)
            w        = spm_superres_lib('get_nii',Nii_w{c});
            Dy       = spm_superres_lib('get_nii',Nii_Dy{c});
            tmp      = sqrt(sum((Dy + w/rho).^2,4));
            z        = (max(tmp - 1/rho,0)./(tmp + eps)).*(Dy + w/rho);
            
            % Update nii
            Nii_z{c} = spm_superres_lib('put_nii',Nii_z{c},z);
        end
        Dy  = [];
        tmp = [];
        w   = [];
        z   = [];
    end        

    %---------------------------
    % Update w
    %---------------------------
    
    parfor (c=1:C,NumWorkers)  
        Dy = spm_superres_lib('get_nii',Nii_Dy{c});
        z  = spm_superres_lib('get_nii',Nii_z{c});
        w  = spm_superres_lib('get_nii',Nii_w{c});
        w  = w + rho*(Dy - z);
        
        % Update nii
        Nii_w{c} = spm_superres_lib('put_nii',Nii_w{c},w);
    end
    mtv = [];
    w   = [];
    z   = [];
    Dy  = [];
    
    %---------------------------
    % Objective function
    %---------------------------
        
    [~,~,dll] = spm_superres_lib('get_ll',Nii_x,Nii_y,Nii_Dy,tau,dat,NumWorkers,Inplane1mm);
    ll        = [ll, dll]; 
%     diff1     = abs((ll(end) - ll(end - 1)));
    diff1     = 2*(ll(end - 1) - ll(end))/(ll(end - 1) + ll(end));
    if Verbose
        fprintf('%3i %10.2f %10.6f %10.6f\n', it, ll(end), diff1, Tolerance);           
    end
    spm_superres_lib('show_stuff',ll,'ll',1,Verbose,ShowZoomed);
    
    %---------------------------
    % Check if converged
    %---------------------------
    
    if it > 20 && ~(diff1 > Tolerance)        
        % Finished
        break
    end
end
toc

if WriteTmpNii    
    % Clean-up temp files
    for c=1:C    
        delete(Nii_z{c}.dat.fname);
        delete(Nii_w{c}.dat.fname);
        delete(Nii_Dy{c}.dat.fname);
        delete(Nii_H{c}.dat.fname);
    end
end

%---------------------------
% Make output
%---------------------------

if isa(Nii_y{1},'nifti')            
    oNii = nifti;
    for c=1:C
        oNii(c) = Nii_y{c};
    end
    
    % Change data-type
    for c=1:C
        y = oNii(c).dat();
        f = Nii_y{c}.dat.fname;
        delete(f);
        
        spm_superres_lib('create_nii',f,y,mat,Nii_x{c}(1).private.dat.dtype,'y', ...
                         Nii_x{c}(1).private.dat.offset,Nii_x{c}(1).private.dat.scl_slope, ...
                         Nii_x{c}(1).private.dat.scl_inter);
    end
else
    oNii = Nii_y;
end
%==========================================================================
