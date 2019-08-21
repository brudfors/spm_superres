function Nii_y = spm_superres(pths,opt)
% Multi-channel total variation (MTV) super-resolution of MR data.
%
% Requires that the SPM software is on the MATLAB path.
% SPM is available from: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
%
% FORMAT Nii_y = spm_superres(pths,opt)
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
% opt   - Algorithm options, more info below.
%
% OUTPUT
% Nii_y - nifti object of C denoised MRIs
%__________________________________________________________________________
% The general principles are described in the following paper:
%
%     Brudfors M, Balbastre Y, Nachev P, Ashburner J.
%     MRI Super-Resolution Using Multi-channel Total Variation.
%     In Annual Conference on Medical Image Understanding and Analysis
%     2018 Jul 9 (pp. 217-228). Springer, Cham.
%
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

if nargin < 1, pths = [];     end
if nargin < 2, opt  = struct; end

% Set up boundary conditions that match the gradient operator
spm_field('boundary',1)

%---------------------------
% Options
%---------------------------

% Algorithm parameters

% scaling of regularisation parameter (lambda)
if ~isfield(opt,'LamScl'),      opt.LamScl      = 10;           end  
% scaling of step-size parameter (rho)
if ~isfield(opt,'RhoScl'),      opt.RhoScl      = 0.5;          end  
% Max number of iterations
if ~isfield(opt,'MaxNiter'),    opt.MaxNiter    = 500;          end  
% Newton iterations
if ~isfield(opt,'NiterNewton'), opt.NiterNewton = 1;            end  
% Convergence threshold
if ~isfield(opt,'Tolerance'),   opt.Tolerance   = 1e-4;         end  
% Run either MTV or indepentend TV denoising
if ~isfield(opt,'DoMTV'),       opt.DoMTV       = true;         end  
% Directory where to write output (and temporary files, which are deleted at end of algorithm)
if ~isfield(opt,'DirOut'),      opt.DirOut      = 'OutputData'; end  
% Clean reference image
if ~isfield(opt,'Nii_y0'),      opt.Nii_y0      = [];           end  
% Number of parfor workers
if ~isfield(opt,'NumWorkers'),  opt.NumWorkers  = Inf;          end  
% Show stuff
if ~isfield(opt,'Verbose'),     opt.Verbose     = true;         end  
% Do preprocessing (register + reslice)
if ~isfield(opt,'DoPreproc'),   opt.DoPreproc   = false;        end  
% Do preprocessing (register)
if ~isfield(opt,'DoCoReg'),     opt.DoCoReg     = false;        end  
% Show one image, zoomed in
if ~isfield(opt,'ShowZoomed'),  opt.ShowZoomed  = true;         end  
% Different test-cases: 0. No testing, 1. brainweb, 2. lena, 3. qmri
if ~isfield(opt,'testcase'),    opt.testcase    = 1;            end  

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
DoPreproc   = opt.DoPreproc;
DoCoReg     = opt.DoCoReg;
ShowZoomed  = opt.ShowZoomed;

%---------------------------
% Get input data
%---------------------------

if opt.testcase && isempty(pths)
    
    %---------------------------
    % Different test-cases (for debugging/testing)
    %---------------------------
        
    imset = struct;
    if     opt.testcase == 1
        imset.ImName = 'brainweb';        
    elseif opt.testcase == 2
        imset.ImName = 'lena';
    elseif opt.testcase == 3
        imset.ImName = 'qmri';
    end   
    imset.LenaNoiseStd = 5*1e1;  % For 'lena', std of additive Gaussian noise
    imset.BrainWebN    = 1;      % For 'brainweb', number of observations of each channel
    imset.BrainWeb2D   = true;  % For 'brainweb', use 2D data
    imset.qmriNumRuns  = 5;      % For 'qmri', number of runs to use
    imset.qmriNumC     = 22;     % For 'qmri', number of channels to use
    imset.qmri2D       = true;   % For 'qmri', use 2D data
    DoPreproc          = false;
    
    DirOut = fullfile(DirOut,imset.ImName);
    if ~(exist(DirOut,'dir') == 7) 
        mkdir(DirOut);
    end

    [Nii_x,Nii_y0] = load_testdata(imset);
else
    %---------------------------
    % Real data
    %---------------------------
    
    % Parse input
    Nii_x = parse_input(pths);
            
    if DoPreproc
        % Register and reslice input
        Nii_x = preproc_input(Nii_x,DirOut,NumWorkers);    
    end
end

%---------------------------
% Estimate rigid alignment matrices
%---------------------------

R = coreg_input(Nii_x,DoCoReg);
    
%---------------------------
% Get projection matrix struct and image properties
%---------------------------

[dat,dm,mat,vx] = get_dat(Nii_x,R);
    
C  = numel(Nii_x); % Number of channels
nm = prod(dm);     % Number of voxels     
N0 = 0;            % Total number of observations
for c=1:C
    N0 = N0 + numel(Nii_x{c});
end

% Show input
show_stuff(Nii_y0,'y',1,Verbose,ShowZoomed);
if dat(1).A(1).do_pm
    show_stuff(Nii_x, ['x (C=' num2str(C) ' N=' num2str(N0) ')'],2,Verbose,true);
else
    show_stuff(Nii_x, ['x (C=' num2str(C) ' N=' num2str(N0) ')'],2,Verbose,ShowZoomed);
end

%---------------------------
% Estimate model parameters
%---------------------------

[tau,lam,rho] = estimate_model_parameters(Nii_x,LamScl,RhoScl,NumWorkers,Verbose);

%---------------------------
% Init variables (creates a load of niftis that will be deleted at end of
% algorithm)
%---------------------------
    
[Nii_y,Nii_z,Nii_w,Nii_Dy,Nii_H] = alloc_vars(Nii_x,dm,mat,DirOut,Verbose);

%---------------------------
% Run algorithm 
%---------------------------

if Verbose
    fprintf('\nRunning %d iterations:\n', MaxNiter);
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
        for it1=1:NiterNewton
                    
            y = get_nii(Nii_y{c});
            z = get_nii(Nii_z{c});
            w = get_nii(Nii_w{c});
           
            % Compute gradient and Hessian
            gr = single(0);
            H  = single(0);
            
            % Conditional part (mask)
            for n=1:N                
                % Mask missing data
                x     = get_nii(Nii_x{c}(n));
                msk   = get_msk(x);
                datcn = dat(c).A(n);
                
                % Gradient part
                gr1       = tau{c}(n)*pm('At',(pm('A',y,datcn) - x),datcn);
                gr1(~msk) = 0;
                gr        = gr + gr1;
                x         = [];
                
                % Hessian part
                if it == 1 && NiterNewton == 1
                    % Compute only once (does not change)
                    H1       = tau{c}(n)*pm('At',pm('A',ones(dm,'single'),datcn),datcn);
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
                Nii_H{c} = put_nii(Nii_H{c},H);
            else
                % Load Hessian
                H = get_nii(Nii_H{c});
            end
            
            gr = gr + diffoperator(w - rho*z,dm,vx,lam{c},'Dt');
            z  = [];
            w  = [];
            
            gr = gr + spm_field('vel2mom',y,[vx 0 rho*lam{c}^2 0]);

            % Do update with spm_field
            y = y - spm_field(H,gr,[vx 0 rho*lam{c}^2 0  1 1]);
            g = []; 
            H = [];
            
            % Update nii
            Nii_y{c} = put_nii(Nii_y{c},y);
        end

        % Compute Dy
        Dy        = diffoperator(y,dm,vx,lam{c},'D');
        Nii_Dy{c} = put_nii(Nii_Dy{c},Dy);
        y         = [];
        Dy        = [];
    end    
    
    show_stuff(Nii_y,'yhat',4,Verbose,ShowZoomed);
    
    %---------------------------
    % Update z
    %---------------------------
    
    if DoMTV
        % MTV
        znorm = single(0);
        parfor (c=1:C,NumWorkers)     
            w     = get_nii(Nii_w{c});
            Dy    = get_nii(Nii_Dy{c});
            znorm = znorm + sum((Dy + w/rho).^2,4);        
        end
        Dy = [];
        w  = [];
        
        znorm = sqrt(znorm);
        mtv   = max(znorm - 1/rho,0)./(znorm + eps);
        znorm = [];

        parfor (c=1:C,NumWorkers)    
            w  = get_nii(Nii_w{c});
            Dy = get_nii(Nii_Dy{c});
            z  = mtv.*(Dy + w/rho);
            
            % Update nii
            Nii_z{c} = put_nii(Nii_z{c},z);
        end    
        Dy = [];
        w  = [];
        z  = [];
        
        show_stuff(mtv,'mtv',3,Verbose,ShowZoomed);
    else
        % Regular TV
        parfor (c=1:C,NumWorkers)
            w        = get_nii(Nii_w{c});
            Dy       = get_nii(Nii_Dy{c});
            tmp      = sqrt(sum((Dy + w/rho).^2,4));
            z        = (max(tmp - 1/rho,0)./(tmp + eps)).*(Dy + w/rho);
            
            % Update nii
            Nii_z{c} = put_nii(Nii_z{c},z);
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
        Dy = get_nii(Nii_Dy{c});
        z  = get_nii(Nii_z{c});
        w  = get_nii(Nii_w{c});
        w  = w + rho*(Dy - z);
        
        % Update nii
        Nii_w{c} = put_nii(Nii_w{c},w);
    end
    mtv = [];
    w   = [];
    z   = [];
    Dy  = [];
    
    %---------------------------
    % Objective function
    %---------------------------
        
    [~,~,dll] = get_ll(Nii_x,Nii_y,Nii_Dy,tau,dat,NumWorkers);
    ll        = [ll, dll]; 
    diff1     = abs((ll(end) - ll(end - 1)));
    if Verbose
        fprintf('%3i %10.2f %10.2f %10.2f\n', it, ll(end), diff1, nm*Tolerance);           
    end
    show_stuff(ll,'ll',1,Verbose,ShowZoomed);
    
    %---------------------------
    % Check if converged
    %---------------------------
    
    if it > 20 && ~(diff1 > Tolerance*nm)        
        % Finished
        break
    end
end
toc

if ~isempty(Nii_y0)
    % For 'lena', show RGB images
    show_stuff(Nii_y0,'y',1,Verbose,ShowZoomed,true);
    show_stuff(Nii_x, 'x',2,Verbose,ShowZoomed,true);
    show_stuff(Nii_y, 'yhat',4,Verbose,ShowZoomed,true);
end

%---------------------------
% Clean-up
%---------------------------
  
for c=1:C    
    delete(Nii_z{c}.dat.fname);
    delete(Nii_w{c}.dat.fname);
    delete(Nii_Dy{c}.dat.fname);
    delete(Nii_H{c}.dat.fname);
    if DoPreproc
        N = numel(Nii_x{c});
        for n=1:N
            delete(Nii_x{c}(n).dat.fname);
        end
    end
end
%==========================================================================

%==========================================================================
% alloc_vars()
function [Nii_y,Nii_z,Nii_w,Nii_Dy,Nii_H] = alloc_vars(Nii_x,dm,mat,DirOut,Verbose)
if Verbose, fprintf('Allocating niftis...'); end

if ~(exist(DirOut,'dir') == 7) 
    mkdir(DirOut);
end

C = numel(Nii_x);

Nii_y  = {nifti};
Nii_z  = {nifti};
Nii_w  = {nifti};    
Nii_Dy = {nifti};   
Nii_H  = {nifti}; 
for c=1:C
    if isa(Nii_x{c}(1),'nifti')        
        f = Nii_x{c}(1).dat.fname;        
    else
        f = ['x' num2str(c) '.nii'];
    end
    [pth,nam] = fileparts(f);
            
    if ~isempty(DirOut)
        pth = DirOut;
    end
    
    fname_y  = fullfile(pth,['y'  nam '.nii']);
    fname_z  = fullfile(pth,['z'  nam '.nii']); 
    fname_w  = fullfile(pth,['w'  nam '.nii']);
    fname_Dy = fullfile(pth,['Dy' nam '.nii']);
    fname_H  = fullfile(pth,['H'  nam '.nii']);

    create_nii(fname_y,zeros(  dm(1:3),   'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
    create_nii(fname_z,zeros( [dm(1:3) 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'z');
    create_nii(fname_w,zeros( [dm(1:3) 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
    create_nii(fname_Dy,zeros([dm(1:3) 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'Dy');
    create_nii(fname_H,zeros(  dm(1:3),   'single'),mat,[spm_type('float32') spm_platform('bigend')],'H');

    Nii_y{c}  = nifti(fname_y);
    Nii_z{c}  = nifti(fname_z);
    Nii_w{c}  = nifti(fname_w);
    Nii_Dy{c} = nifti(fname_Dy);
    Nii_H{c}  = nifti(fname_H);
end
if Verbose, fprintf('done!\n'); end
%==========================================================================

%==========================================================================
% create_nii()
function Nii = create_nii(pth,dat,mat,dtype,descrip,offset,scl_slope,scl_inter)
if nargin<6, offset    = 0; end
if nargin<7, scl_slope = 1; end
if nargin<8, scl_inter = 0; end

if exist(pth,'file')==2, delete(pth); end

Nii         = nifti;
dm          = size(dat);
Nii.dat     = file_array(pth,dm,dtype,offset,scl_slope,scl_inter);
Nii.mat     = mat;
Nii.mat0    = mat;
Nii.descrip = descrip;
create(Nii);

Nii.dat(:) = dat(:);
%==========================================================================

%==========================================================================
% diffoperator()
function out = diffoperator(in,dm,vx,lam,type)
Dx  = sparse(toeplitz([-1 1 zeros(1,dm(1)-2)],[-1 zeros(1,dm(1)-1)]))/vx(1); Dx(1,1)=0;
Dy  = sparse(toeplitz([-1 1 zeros(1,dm(2)-2)],[-1 zeros(1,dm(2)-1)]))/vx(2); Dy(1,1)=0;
if dm(3) == 1
    Dz = sparse(0);
else
    Dz = sparse(toeplitz([-1 1 zeros(1,dm(3)-2)],[-1 zeros(1,dm(3)-1)]))/vx(3); Dz(1,1)=0;
end

D = lam*[kron(speye(dm(3)),kron(speye(dm(2)),Dx          ))
         kron(speye(dm(3)),kron(Dy,          speye(dm(1))))
         kron(Dz,          kron(speye(dm(2)),speye(dm(1))))];
clear Dx Dy Dz

if strcmp(type,'D')
    out = full(reshape(D*double(in(:)),[dm(1:3) 3]));
elseif strcmp(type,'Dt')
    out = full(reshape(D'*double(in(:)),dm(1:3)));
else
    error('Input error!');
end
out = single(out);
%==========================================================================

%==========================================================================
% coreg_input()
function R = coreg_input(Nii_x,DoCoReg,ref)
if nargin < 3, ref = [1 1]; end

C  = numel(Nii_x);
Vr = spm_vol(Nii_x{ref(1)}(ref(2)).dat.fname);
R  = cell(1,C);
for c=1:C
    N    = numel(Nii_x{c}); 
    R{c} = cell(1,N);
    for n=1:N
        if (c == ref(1) && n == ref(2)) || ~DoCoReg
            R{c}{n} = eye(4);
        else    
            Vm      = spm_vol(Nii_x{c}(n).dat.fname);
            x       = spm_coreg(Vr,Vm);
            R{c}{n} = spm_matrix(x(:)');
        end        
    end
end
%==========================================================================

%==========================================================================
% estimate_model_parameters()
function [tau,lam,rho] = estimate_model_parameters(Nii,LamScl,RhoScl,NumWorkers,Verbose,LenaNoiseStd)
if nargin < 6, LenaNoiseStd = 1e1; end   % For 'lena', std of additive Gaussian noise

C   = numel(Nii);
tau = cell(1,C); 
lam = cell(1,C);
parfor (c=1:C,NumWorkers)
    N = numel(Nii{c});
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

% This value of rho seems to lead to reasonably good convergence
rho = RhoScl*sqrt(mean(reshape(cell2mat(tau),[],1)))/mean(cell2mat(lam));
if Verbose
    fprintf('step-size -> rho=%f\n', rho);
end
%==========================================================================

%==========================================================================
% get_dat()
function [dat,dm,mat,vx] = get_dat(Nii_x,R)

%---------------------------
% Check if projection matrices need to be used
%---------------------------

C  = numel(Nii_x);
vx = [];
dm = [];
for c=1:C
    N = numel(Nii_x{c});
    for n=1:N
        vxc = sqrt(sum(Nii_x{c}(n).mat(1:3,1:3).^2)); % Voxel sizes   
        dmc = size(Nii_x{c}(n).dat);                  % Image size
        dmc = [dmc 1];

        vx = [vx; round(vxc,2)];
        dm = [dm; dmc];
    end
end

vx = unique(vx,'rows');
dm = unique(dm,'rows');

%---------------------------
% Build projection matrix struct (dat)
%---------------------------

if size(dm,1) == 1 && all(vx == 1)
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
    [mat,dm,vx] = max_bb_orient(Nii_x,R);    
    dat         = init_dat(Nii_x,mat,dm,R);
end
%==========================================================================

%==========================================================================
% get_ll()
function [ll1,ll2,ll] = get_ll(Nii_x,Nii_y,Nii_Dy,tau,dat,NumWorkers)
% log posterior
C   = numel(Nii_y);
ll1 = 0;
ll2 = 0;
parfor (c=1:C,NumWorkers)
    y = get_nii(Nii_y{c});
    for n=1:numel(Nii_x{c})
        x   = get_nii(Nii_x{c}(n));
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
%==========================================================================

%==========================================================================
% get_msk()
function msk = get_msk(in)
msk = isfinite(in) & in ~= 0;
%==========================================================================

%==========================================================================
% get_nii()
function img = get_nii(nii)
img = single(nii.dat());
%==========================================================================

%==========================================================================
% load_testdata()
function [Nii_x,Nii_y0] = load_testdata(imset)

ImName       = imset.ImName;
LenaNoiseStd = imset.LenaNoiseStd;
BrainWebN    = imset.BrainWebN;
BrainWeb2D   = imset.BrainWeb2D;
qmriNumRuns  = imset.qmriNumRuns;
qmriNumC     = imset.qmriNumC;
qmri2D       = imset.qmri2D;

if strcmp(ImName,'lena')
    
    %----------------
    % Load Lena image (y), add some noise to it plus remove some parts (x)
    %----------------
    
    % Repeatable random numbers..
    rng('default');
    rng(1);

    
    
    im = imread('TestData/lena/lena_std.tif');
    im = single(im);
    C  = size(im,3);
        
    y  = cell(1,C);
    for c=1:C
        y{c} = im(:,:,c);%/max(im(:));
    end

    % Add noise
    x = cell(1,C);
    n = 1;
    for c=1:C
        x{c}{n}              = y{c} + LenaNoiseStd*randn(size(y{c}));
        x{c}{n}(x{c}{n} < 0) = 0;
    end                             % Number of voxels
    
    % Missing data
    dm = size(x{c}{n});
    dm = [dm 1];
    x{1}{n}(:,1:ceil(dm(2)/3))                = NaN;
    x{2}{n}(:,floor(dm(2)/3):2*ceil(dm(2)/3)) = NaN;
    x{3}{n}(:,2*floor(dm(2)/3):end)           = NaN;    
    
    Nii_x  = cell(1,C);
    Nii_y0 = cell(1,C);
    for c=1:C        
        Nii_x{c}(1).mat = eye(4);
        Nii_x{c}(n).dat = x{c}{n};
        Nii_y0{c}.dat   = y{c};
    end         
    
elseif strcmp(ImName,'brainweb')
    
    %----------------
    % BrainWeb multi-channel
    %----------------
    
    DirRef = 'TestData/brainweb';    
    C      = 3;
    
    Nii_x = cell(1,C);
    for n=1:BrainWebN
        DirSim        = fullfile('./SimulatedData/brainweb',['n' num2str(n)]);
        [Nii3d,Nii2d] = SimulateData('DirRef',DirRef, 'DirSim',DirSim, ...
                                     'Random',[true false false false false true]);
        for c=1:C
            if n == 1
                Nii_x{c} = nifti;
            end
            if BrainWeb2D
                Nii_x{c}(n) = Nii2d(c);
            else
                Nii_x{c}(n) = Nii3d(c);
            end
        end
    end 
    Nii_y0 = [];
elseif strcmp(ImName,'qmri')    
    
    %----------------
    % qMRI
    %----------------
    
    if qmri2D
        load('TestData/qmri_reg_2d.mat');
        Nii_x = Nii_x(1:qmriNumC);
        for c=1:qmriNumC  
            Nii_x{c} = Nii_x{c}(1:qmriNumRuns);            
        end
    else
        DirRef = 'TestData/qmri-reg';                          
        DirSim = 'TestData/qmri-2D';    

        Files   = cell(1,qmriNumRuns);
        if exist(DirSim,'dir') == 7, rmdir(DirSim,'s'); end; mkdir(DirSim);
        for i=1:qmriNumRuns
            dir_run  = fullfile(DirRef,['run_0' num2str(i)]);        

            if qmri2D
                Dir2Dr  = fullfile(DirSim,['run_0' num2str(i)]);
                dir_run = make2d(dir_run,Dir2Dr,2);
            end

            Files{i} = spm_select('FPList',dir_run,'^.*\.nii$');
        end

        C     = min(qmriNumC,size(Files{1},1));
        Nii_x = cell(1,C);    
        for c=1:C
            Nii_x{c} = nifti;    
            for n=1:qmriNumRuns
                Nii_x{c}(n) = nifti(strtrim(Files{n}(c,:)));
            end    
        end
    end
    Nii_y0 = [];
else
    error('Unknown type!')
end
%==========================================================================

%==========================================================================
% max_bb_orient()
function [mat,dm,vx] = max_bb_orient(Nii,R,vx)

if nargin < 3, vx = [1 1 1]; end

mn = [ Inf  Inf  Inf]';
mx = [-Inf -Inf -Inf]';
for c=1:numel(Nii)
    N = numel(Nii{c});
    for n=1:N
        dmn = size(Nii{c}(n).dat);
        
        if numel(dmn) == 2
            dmn(3) = 0; 
        end

        t = uint8(0:7);
        s = diag(dmn+1)*double([bitshift(bitand(t,bitshift(uint8(1),1-1)),1-1)
                              bitshift(bitand(t,bitshift(uint8(1),2-1)),1-2)
                              bitshift(bitand(t,bitshift(uint8(1),3-1)),1-3)]);
                          
        mat = R{c}{n}\Nii{c}(n).mat;
        
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
%==========================================================================

%==========================================================================
% parse_input()
function Nii = parse_input(pths)
C   = numel(pths);
Nii = cell(1,C);
for c=1:C
    Nii{c} = nifti(pths{c});
end
%==========================================================================

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
%==========================================================================

%==========================================================================
% pm()
function varargout = pm(nam,varargin)
if strcmp(nam,'A')
    %---------------------------
    % Forward operator
    %---------------------------

    % Parse input
    y = varargin{1};
    A = varargin{2};
    
    if A.do_pm
        x = apply_proj(y,A,nam);
    else
        % Identity
        x = y;
    end
    
    % Make output
    varargout{1} = x;
elseif strcmp(nam,'At')
    %---------------------------
    % Adjoint operator
    %---------------------------
    
    % Parse input
    x = varargin{1};
    A = varargin{2};
    
    if A.do_pm
        y = apply_proj(x,A,nam);
    else
        % Identity
        y = x;
    end
    
    % Make output
    varargout{1} = y;
else
    error('Input argument error!')
end
%==========================================================================

%==========================================================================
% put_nii()
function nii = put_nii(nii,img)
nii.dat(:) = img(:);
%==========================================================================

%==========================================================================
% show_stuff()
function show_stuff(in,nam,np,verbose,ShowZoomed,isrgb,nr,nc,figname)
if nargin < 4, verbose    = true; end
if nargin < 5, ShowZoomed = false; end
if nargin < 6, isrgb      = false; end
if nargin < 7, nr         = 2; end
if nargin < 8, nc         = 2; end
if nargin < 9, figname    = '(SPM) spm_denoise'; end

fig                  = findobj('Type', 'Figure', 'Name', figname);
if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
set(0, 'CurrentFigure', fig);      
    
if isempty(in)
    clf(fig);
    return
end

if verbose  
    subplot(nr,nc,np);
    if strcmp(nam,'ll')        
        plot(in,'LineWidth',2);    
    else
        if ShowZoomed
            if iscell(in)
                im = get_nii(in{1}(1));
                mx = max(im(:));                
                if mx == 0
                    mx = 1;
                end
                im = im/mx;
            else
                mx = max(in(:));
                if mx == 0
                    mx = 1;
                end
                im = in/mx;
            end
            
            dm = size(im);
            val = 0.25;
            y  = floor(val*dm(2):dm(2) - val*dm(2));
            x  = floor(val*dm(1):dm(1) - val*dm(1));
            im = im(x,y,:);
        else
            if iscell(in)
                C  = numel(in);
                im = [];
                for c=1:C                
                    N = numel(in{c});
                    for n=1:N                
                        imn = get_nii(in{c}(n));
                        mx  = max(imn(:));                
                        if mx == 0
                            mx = 1;
                        end
                        im = cat(4,im,imn/mx);                    
                    end
                end                
            else
                mx = max(in(:));
                if mx == 0
                    mx = 1;
                end
                im = in/mx;
            end
        end
        
        dm = size(im);
        dm = [dm 1];
        z  = round(dm(3)/2);
        if isrgb
            imshow(squeeze(im(:,:,z,:)));
        else
            if size(im,3) == 1 && size(im,4) == 1
                imagesc(im);
            else
                montage(squeeze(im(:,:,z,:)),'DisplayRange',[0 max(im(:)) + eps]);
            end
            colormap(gray);
        end
        axis off image;
    end
    title(nam)
    drawnow
end
%==========================================================================