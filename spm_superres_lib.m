function varargout = spm_superres_lib(varargin)
%__________________________________________________________________________
% Library of functions for spm_superres.
%
% FORMAT out = spm_superres_lib('name'in)
%
% FORMAT help spm_superres_lib>function
% Returns the help file of the selected function.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

if nargin == 0
    help spm_superres_lib
    error('Not enough argument. Type ''help spm_superres_lib'' for help.');
end
id = varargin{1};
varargin = varargin(2:end);
switch lower(id) 
    case 'alloc_vars'
        [varargout{1:nargout}] = alloc_vars(varargin{:});     
    case 'apply_affine'
        [varargout{1:nargout}] = apply_affine(varargin{:});                           
    case 'check_do_write'
        [varargout{1:nargout}] = check_do_write(varargin{:});          
    case 'coreg_input'
        [varargout{1:nargout}] = coreg_input(varargin{:});                  
    case 'create_nii'
        [varargout{1:nargout}] = create_nii(varargin{:});                      
    case 'diffoperator'
        [varargout{1:nargout}] = diffoperator(varargin{:});                              
    case 'estimate_model_parameters'
        [varargout{1:nargout}] = estimate_model_parameters(varargin{:});                                      
    case 'get_dat'
        [varargout{1:nargout}] = get_dat(varargin{:});              
    case 'get_ll'
        [varargout{1:nargout}] = get_ll(varargin{:});              
    case 'get_msk'
        [varargout{1:nargout}] = get_msk(varargin{:});        
    case 'get_nii'
        [varargout{1:nargout}] = get_nii(varargin{:});          
    case 'get_opt'
            [varargout{1:nargout}] = get_opt(varargin{:});             
    case 'parse_input'
            [varargout{1:nargout}] = parse_input(varargin{:});          
    case 'pm'
            [varargout{1:nargout}] = pm(varargin{:});                      
    case 'put_nii'
            [varargout{1:nargout}] = put_nii(varargin{:});                          
    case 'show_stuff'
            [varargout{1:nargout}] = show_stuff(varargin{:});                          
    otherwise
        help spm_superres_lib
        error('Unknown function %s. Type ''help spm_superres_lib'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% alloc_vars()
function [Nii_y,Nii_z,Nii_w,Nii_Dy,Nii_H] = alloc_vars(WriteTmpNii,Nii_x,dm,mat,DirOut,Verbose)
if Verbose, fprintf('Allocating niftis...'); end

if ~(exist(DirOut,'dir') == 7) && ~isempty(DirOut)
    mkdir(DirOut);
end

C = numel(Nii_x);

Nii_y  = {nifti};
if WriteTmpNii
    Nii_z  = {nifti};
    Nii_w  = {nifti};    
    Nii_Dy = {nifti};   
    Nii_H  = {nifti}; 
else
    Nii_z  = {struct};
    Nii_w  = {struct};    
    Nii_Dy = {struct};   
    Nii_H  = {struct}; 
end
for c=1:C
    if isstruct(Nii_x{c}(1))        
        f = Nii_x{c}(1).fname;        
    else
        f = ['img' num2str(c) '.nii'];
    end
    [pth,nam] = fileparts(char(f));
            
    if ~isempty(DirOut)
        pth = DirOut;
    end
    
    fname_y  = fullfile(pth,['y'  nam '.nii']);
    create_nii(fname_y,zeros(dm(1:3),'single'),mat,[spm_type('float32') spm_platform('bigend')],'y');
    Nii_y{c} = nifti(fname_y);

    if WriteTmpNii
        fname_z  = fullfile(pth,['z'  nam '.nii']); 
        fname_w  = fullfile(pth,['w'  nam '.nii']);
        fname_Dy = fullfile(pth,['Dy' nam '.nii']);
        fname_H  = fullfile(pth,['H'  nam '.nii']);

        create_nii(fname_z,zeros( [dm(1:3) 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'z');
        create_nii(fname_w,zeros( [dm(1:3) 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'w');
        create_nii(fname_Dy,zeros([dm(1:3) 3],'single'),mat,[spm_type('float32') spm_platform('bigend')],'Dy');
        create_nii(fname_H,zeros(  dm(1:3),   'single'),mat,[spm_type('float32') spm_platform('bigend')],'H');

        Nii_z{c}  = nifti(fname_z);
        Nii_w{c}  = nifti(fname_w);
        Nii_Dy{c} = nifti(fname_Dy);
        Nii_H{c}  = nifti(fname_H);
    else
        Nii_z{c}.dat  = zeros([dm(1:3) 3],'single');
        Nii_w{c}.dat  = zeros([dm(1:3) 3],'single');
        Nii_Dy{c}.dat = zeros([dm(1:3) 3],'single');
        Nii_H{c}.dat  = zeros(dm(1:3),'single');
    end
end
if Verbose, fprintf('done!\n'); end
end
%==========================================================================

%==========================================================================
% apply_affine()
function y = apply_affine(M,dm)
[x0,y0,z0] = ndgrid(single(1:dm(1)),...
                    single(1:dm(2)),...
                    single(1:dm(3)));
y          = cat(4,M(1,1)*x0 + M(1,2)*y0 + M(1,3)*z0 + M(1,4), ...
                   M(2,1)*x0 + M(2,2)*y0 + M(2,3)*z0 + M(2,4), ...
                   M(3,1)*x0 + M(3,2)*y0 + M(3,3)*z0 + M(3,4));
if dm(3) == 1
    y(:,:,:,end) = 1;
end
end
%========================================================================== 

%==========================================================================
% check_do_write()
function do_write = check_do_write(C,dm,memmx,nw)
% z  C*dm*3*4b
% w  C*dm*3*4b
% Dy C*dm*3*4b
% H  C*dm*4b
% and X and Y
nm     = prod(dm);
unt    = 4; % One float is four bytes
memreq = C*nm*3*unt + C*nm*3*unt + C*nm*3*unt + C*nm*unt + 2*min(min(nw,32),C)*nm*unt;
memreq = memreq/1e6;
if memreq > memmx
    do_write = true;
else
    do_write = false;
end
end
%==========================================================================

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

if Nii_x{1}(1).dim(3) > 1
    
    Vr = Nii_x{ref(1)}(ref(2));
    for c=1:C
        N = numel(Nii_x{c});         
        for n=1:N
            if (c == ref(1) && n == ref(2)) || ~DoCoReg
                continue
            else    
                Vm      = Nii_x{c}(n);
                x       = spm_coreg(Vr,Vm);
                R{c}{n} = spm_matrix(x(:)');
            end        
        end
    end
end
end
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
end
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
end
%==========================================================================

%==========================================================================
% estimate_model_parameters()
function [tau,lam,rho] = estimate_model_parameters(Nii,LamScl,RhoScl,NumWorkers,Verbose)

C    = numel(Nii);
tau  = cell(1,C); 
lam  = cell(1,C);
info = cell(1,C);
parfor (c=1:C,NumWorkers)
            
    % Estimate MR image noise
    [NoiseStd,mu,info{c}] = noise_estimate(Nii{c});
    tau{c}                = 1./(NoiseStd.^2);                                  
    
    lam{c} = LamScl/double(mean(mu));
    
%     % Scale with number of observations to give a little more weight to the
%     % prior when there are multiple observations of the same channel
%     lam{c} = sqrt(N)*lam{c}; 
    
    if Verbose
        N = numel(Nii{c});
        for n=1:N
            fprintf('c=%i, i=%i | sd=%f, mu=%f -> tau=%f, lam=%f\n', c, n, NoiseStd(n), mu(n), tau{c}(n), lam{c});
        end
    end
end

if Verbose >= 2
    
    %---------------------------
    % Show fit(s)
    %---------------------------
    
    figname              = '(spm_superres) Parameter estimates';
    fig                  = findobj('Type', 'Figure', 'Name', figname);
    if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
    set(0, 'CurrentFigure', fig);      
    
    N0 = 0;
    for c=1:C
        N = numel(Nii{c});
        for n=1:N
            N0 = N0 + 1;
        end
    end
    nr = floor(sqrt(N0));
    nc = ceil(N0/nr);  
    
    cnt = 1;
    for c=1:C       
        N = numel(Nii{c});
        for n=1:N
            subplot(nr,nc,cnt);                        
            plot(info{c}(n).x(:),info{c}(n).h/sum(info{c}(n).h)/info{c}(n).md,'b.', ...
                 info{c}(n).x(:),info{c}(n).sp,'r', ...
                 info{c}(n).x(:),info{c}(n).p,'--'); 
             
%             FontSize = 18;
%             set(gcf, 'Color', 'w')
%             set(gca,'FontSize',FontSize)
%             legend('Empirical','Mixture fit','Air class','Brain class')
%             xlabel('Image intensity')
%             title('RMM Fitted to MRI Intensity Histogram')
    
            cnt = cnt + 1;
        end
    end
    drawnow
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

        vx = [vx; round(vxc*10^2)/10^2];
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

%==========================================================================
% get_msk()
function msk = get_msk(in)
msk = isfinite(in) & in ~= 0;
end
%==========================================================================

%==========================================================================
% get_nii()
function img = get_nii(nii)
if isa(nii, 'nifti') || (isstruct(nii) && ~isfield(nii, 'private'))
      nii = nii.dat;
end
if isstruct(nii)
      img = single(spm_read_vols(nii));
else
      [idx{1:ndims(nii)}] = deal(':');
      img = single(subsref(nii,struct('type','()','subs',{idx})));
end
% if isa(nii, 'nifti')
%     img = single(nii.dat());
% elseif isstruct(nii)
%     if isfield(nii, 'private')
%         img = single(spm_read_vols(nii));
%     else
%         img = single(nii.dat());
%     end
% else
%     img = single(nii());
% end
end
%==========================================================================

%==========================================================================
% get_opt()
function opt = get_opt(opt)
% scaling of regularisation parameter (lambda)
if ~isfield(opt,'LamScl'),      opt.LamScl      = 10;           end  
% scaling of step-size parameter (rho)
if ~isfield(opt,'RhoScl'),      opt.RhoScl      = 1;            end  
% Max number of iterations
if ~isfield(opt,'MaxNiter'),    opt.MaxNiter    = 30;           end  
% Newton iterations
if ~isfield(opt,'NiterNewton'), opt.NiterNewton = 1;            end  
% Convergence threshold
if ~isfield(opt,'Tolerance'),   opt.Tolerance   = 1e-4;         end  
% Run either MTV or indepentend TV denoising
if ~isfield(opt,'DoMTV'),       opt.DoMTV       = true;         end  
% Directory where to write output (and temporary files, which are deleted at end of algorithm)
if ~isfield(opt,'DirOut'),      opt.DirOut      = ''; end  
% Clean reference image
if ~isfield(opt,'Nii_y0'),      opt.Nii_y0      = [];           end  
% Number of parfor workers
if ~isfield(opt,'NumWorkers'),  opt.NumWorkers  = 8;            end  
% Show stuff
if ~isfield(opt,'Verbose'),     opt.Verbose     = 1;            end  
% Do preprocessing (register)
if ~isfield(opt,'DoCoReg'),     opt.DoCoReg     = true;         end  
% Show one image, zoomed in
if ~isfield(opt,'ShowZoomed'),  opt.ShowZoomed  = false;        end  
% Memory limit for when to allocate variables as niftis
if ~isfield(opt,'MaxMem'),      opt.MaxMem      = 4096;         end 
% Reconstruction voxel size, if 0, set to smallest available
if ~isfield(opt,'VoxSize'),     opt.VoxSize     = 1;            end 
% Downsample inplane resolution to 1 mm
if ~isfield(opt,'Inplane1mm'),  opt.Inplane1mm  = true;         end 
% Do just denoising, without super-resolving
if ~isfield(opt,'Denoise'),     opt.Denoise     = false;        end 
end
%==========================================================================

%==========================================================================
% parse_input()
function Nii = parse_input(pths)
if iscell(pths)
    C   = numel(pths);
    Nii = cell(1,C);
    for c=1:C
        Nii{c} = spm_vol(pths{c});
    end
elseif ischar(pths) && (exist(pths,'dir') == 7)
    dirs = spm_select('List',pths,'dir');    
    C    = size(dirs,1);
    if C > 0
        Nii = cell(1,C);
        for c=1:C
            d      = deblank(dirs(c,:));
            d      = fullfile(pths,d);
            files  = spm_select('FPList',d,'^.*\.nii$');
            Nii{c} = spm_vol(files);
        end
    else
        files = spm_select('FPList',pths,'^.*\.nii$');
        C     = size(files,1);
        Nii   = cell(1,C);
        for c=1:C
            f      = deblank(files(c,:));
            Nii{c} = spm_vol(f);
        end
    end
else
    error('Input error!')
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
end
%==========================================================================

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
end
%==========================================================================

%==========================================================================
% show_stuff()
function show_stuff(in,nam,np,verbose,ShowZoomed,isrgb,nr,nc,figname)
if nargin < 4, verbose    = 2; end
if nargin < 5, ShowZoomed = false; end
if nargin < 6, isrgb      = false; end
if nargin < 7, nr         = 2; end
if nargin < 8, nc         = 2; end
if nargin < 9, figname    = '(spm_superres) Algorithm progress'; end

if verbose >= 2
    fig                  = findobj('Type', 'Figure', 'Name', figname);
    if isempty(fig), fig = figure('Name', figname, 'NumberTitle', 'off'); end
    set(0, 'CurrentFigure', fig);      

    if isempty(in)
        clf(fig);
        return
    end

    subplot(nr,nc,np);
    if strcmp(nam,'ll')        
        plot(in(min(3,numel(in)):end),'LineWidth',2);    
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
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% apply_proj()
function out = apply_proj(in,A,method)

spm_diffeo('bound',1); % OK?
spm_field('bound',1);  % match the gradient operator

% Projection parameters
R    = A.R;
Mmu  = A.mat0;
Mf   = A.mat;
dmmu = A.dm0;
dmf  = A.dm;                
vsmu = sqrt(sum(Mmu(1:3,1:3).^2));        
vsf  = sqrt(sum(Mf(1:3,1:3).^2));
scl  = abs(det(Mmu(1:3,1:3)))^(1/3);

samp           = vsf./scl;    
samp(samp < 1) = 1;
D              = diag([samp 1]);        
% smo            = sqrt(max(scl.^2-vsf.^2,0))*sqrt(8*log(2));
smo            = vsf./scl; 
smo(smo <= 1)  = 0;
[~,ix]         = max(smo);
gap            = 1/3;
% smo(ix)        = smo(ix) - gap*smo(ix);
sd             = smo./(2*sqrt(2*log(2)));
sdscl          = 4;

off         = -(sdscl*sd)';
Moff        = eye(4);
Moff(1:3,4) = off;
dmoff       = ceil((D(1:3,1:3)*dmf(1:3)')' + 2*sdscl.*sd);
% dmoff       = ceil((D(1:3,1:3)*dmf(1:3)')');

% Verbose options
verbose = false;
nr      = 2;
nc      = 2;
plotnum = 1;
fignum  = 666;

% Do projection, either A or At
if strcmp(method,'A')            
    
    if verbose
        % Verbose
        show_res(in,fignum,nr,nc,plotnum,'in (A)');
        plotnum = plotnum + 1;
    end
    
    % Pull to subject space (with template space resolution)   
    Ty  = Mmu\(R\(Mf/D));    
    Ty  = Ty*Moff;    
    dmy = dmoff;
    out = spm_diffeo('pull',in,apply_affine(Ty,dmy)); 

    if verbose
        % Verbose
        show_res(out,fignum,nr,nc,plotnum,'pull1 (A)');
        plotnum = plotnum + 1;
    end
    
    % Apply slice-profile        
    spm_smooth(out,out,smo);
    
    if verbose
        % Verbose
        show_res(out,fignum,nr,nc,plotnum,'smooth (A)');
        plotnum = plotnum + 1;
    end
    
    % Pull to subject space resolution
    Ty  = D;
    Ty  = Moff\Ty; 
    dmy = dmf;
    out = spm_diffeo('pull',out,apply_affine(Ty,dmy)); 
%     out = out(1:D(1,1):end,1:D(2,2):end,1:D(3,3):end);
        
    if verbose
        % Verbose
        show_res(out,fignum,nr,nc,plotnum,'pull2 (A)');
        plotnum = plotnum + 1;
    end
elseif strcmp(method,'At')
    fignum = fignum + 1;
    
    if verbose
        % Verbose
        show_res(in,fignum,nr,nc,plotnum,'in (At)');
        plotnum = plotnum + 1;
    end
    
    % Push to template space resolution
    Ty  = D;
    Ty  = Moff\Ty; 
    dmy = dmoff;
    out = spm_diffeo('push',in,apply_affine(Ty,dmf),dmy); 
%     out                                     = zeros(dmy,'single');
%     out(1:D(1,1):end,1:D(2,2):end,1:D(3,3)) = in;
    
    if verbose
        % Verbose
        show_res(out,fignum,nr,nc,plotnum,'push1 (At)');
        plotnum = plotnum + 1;
    end
    
    % Apply slice-profile    
    spm_smooth(out,out,smo);

    if verbose
        % Verbose
        show_res(out,fignum,nr,nc,plotnum,'smooth (At)');
        plotnum = plotnum + 1;
    end
    
    % Push to template space (with subject space resolution)
    Ty  = Mmu\(R\(Mf/D));
    Ty  = Ty*Moff;
    dmy = dmoff;
    out = spm_diffeo('push',out,apply_affine(Ty,dmy),dmmu); 
    
    if verbose
        % Verbose
        show_res(out,fignum,nr,nc,plotnum,'push2 (At)');
        plotnum = plotnum + 1;
    end    
end

%---------------------------
% Nested functions
%---------------------------

function show_res(in,fignum,nr,nc,plotnum,nam,ix)
    if nargin < 6, nam = ''; end
    if nargin < 7, ix  = [1 2 3]; end

    figure(fignum);
    subplot(nr,nc,plotnum);
    imagesc(permute(in,ix)); axis off xy image; colormap(gray);
    title(nam);
end

function check_adjoint(dat,n)
    u   = randn(dat.dm, 'single');
    v   = randn(dat.A(n).dm, 'single');
    Au  = apply_proj(u,dat,n,'A');
    Atv = apply_proj(v,dat,n,'At');

    v_dot_Au  = v(:)' * Au(:);
    Atv_dot_v = Atv(:)' * u(:);

    disp(v_dot_Au - Atv_dot_v)
end

end
%==========================================================================

%==========================================================================
% blur_fun()
function f = blur_fun(dm,mat,vx)
if nargin<1, dm = [64 64]; end
if nargin<2, mat = eye(numel(dm)); end
if nargin<3, vx = ones(1,numel(dm)); end

if any(size(mat)~=numel(dm)) || numel(vx)~=numel(dm), error('Incompatible dimensions.'); end

% Grid in frequency space
r        = cell(1,numel(dm));
for i=1:numel(dm) 
    r{i} = single([0:ceil(dm(i)/2-1) -floor(dm(i)/2):-1]'*pi/dm(i)); 
end
X        = cell(1,numel(dm));
[X{:}]   = ndgrid(r{:});
clear r

% Transform
Y            = cell(size(X));
for i=1:numel(dm)
    Y{i}     = single(0);
    for j=1:numel(dm) 
        Y{i} = Y{i} + mat(i,j)*X{j}; 
    end
end
clear X

% Window function
f     = single(0);
for i=1:numel(dm) 
    f = f + Y{i}.^2; 
end    
f     = ((cos(min(f,pi^2/4)*4/pi) + 1)/2);

% Incorporate voxel size
for i=1:numel(dm)
    tmp                 = sin((vx(i))*Y{i})./(Y{i}.*cos(Y{i}/pi^(1/2)));
    tmp(~isfinite(tmp)) = vx(i);
    f                   = f.*tmp;
end
end
%==========================================================================

%==========================================================================
% get_slice_gap()
function gap = get_slice_gap(gap,Nii_x,gapunit)
% Construct slice-gap
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii_x);

if isscalar(gap)
    % Find thick-slice direction automatically
    G = cell(1,C);
    for c=1:C
        N    = numel(Nii_x{c});
        G{c} = cell(1,N);
        for n=1:N
            G{c}{n} = zeros(1,3);
            
            thickslice_dir          = get_thickslice_dir(Nii_x{c}(n));
            G{c}{n}(thickslice_dir) = gap;            
        end
    end    
    gap = G;
else
    % x-, y-, z-directions given
    gap = padarray(gap, [0 max(0,C-numel(gap))], 'replicate', 'post');
    for c=1:C
        if ~iscell(gap{c})
            gap{c} = {gap{c}};
        end
        gap{c} = padarray(gap{c}, [0 max(0,numel(Nii_x{c})-numel(gap{c}))], 'replicate', 'post');
        for i=1:numel(gap{c})
            if isempty(gap{c}{i})
                gap{c}{i} = 0;
            end
            gap{c}{i} = padarray(gap{c}{i}, [0 max(0,3-numel(gap{c}{i}))], 'replicate', 'post');
        end
    end
end

if strcmp(gapunit,'%')
    % Convert from percentage to mm
    for c=1:C
        N    = numel(Nii_x{c});
        for n=1:N       
            mat = Nii_x{c}(n).mat;
            vx  = sqrt(sum(mat(1:3,1:3).^2));   
            
            gap{c}{n} = vx.*gap{c}{n};
        end
    end    
end
end
%==========================================================================

%==========================================================================
% get_thickslice_dir()
function thickslice_dir = get_thickslice_dir(Nii)
mat                = Nii.mat;
vx                 = sqrt(sum(mat(1:3,1:3).^2));   
[~,thickslice_dir] = max(vx);
end
%==========================================================================

%==========================================================================
% get_window()
function window = get_window(window,Nii_x)
% Construct slice-profile
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

C = numel(Nii_x);

if isempty(window)
    % Find thick-slice direction automatically and set defaults:
    % In-plane: Gaussian, Through-plane: Rectangle
    W = cell(1,C);
    for c=1:C
        N    = numel(Nii_x{c});
        W{c} = cell(1,N);
        for n=1:N
            W{c}{n} = ones(1,3); % All Gaussian
            
%             thickslice_dir          = get_thickslice_dir(Nii_x{c}(n));
%             W{c}{n}(thickslice_dir) = 2; % Rectangle
        end
    end    
    window = W;    
else
    % x-, y-, z-directions given
    if ~iscell(window)
        window = {window};
    end
    window = padarray(window, [0 max(0,C-numel(window))], 'replicate', 'post');
    for c=1:C
        if ~iscell(window{c})
            window{c} = {window{c}};
        end
        window{c} = padarray(window{c}, [0 max(0,numel(Nii_x{c})-numel(window{c}))], 'replicate', 'post');
        for i=1:numel(window{c})
            if isempty(window{c}{i})
                window{c}{i} = 2;
            end
            window{c}{i} = padarray(window{c}{i}, [0 max(0,3-numel(window{c}{i}))], 'replicate', 'post');
        end
    end
end
end
%==========================================================================

%==========================================================================
% init_A()
function dat = init_A(Nii,Mmu,dmmu,R,D,window,gap)
% Initialise projection matrices (stored in dat struct)
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging
N          = numel(Nii); % Number of LR images
vs         = sqrt(sum(Mmu(1:3,1:3).^2));
for n=1:N % Loop over LR images
    Mf  = Nii(n).mat;    
    dmf = Nii(n).dim;
    dmf = [dmf 1];

    dat.A(n).mat0 = Mmu;
    dat.A(n).dm0  = dmmu;   

    dat.A(n).mat = Mf/D{n};    
    dat.A(n).dm  = floor(D{n}(1:3,1:3)*dmf(1:3)')';
    dat.A(n).win = window{n};
    dat.A(n).gap = gap{n};
    
    dat.A(n).do_pm = true;
    dat.A(n).R     = R{n};
    dat.A(n).D     = inv(D{n});
    
    M          = Mmu\Mf;
    
%     R          = (M(1:3,1:3)/diag(sqrt(sum(M(1:3,1:3).^2))))';
%     dat.A(n).S = blur_fun(dm,R,sqrt(sum(M(1:3,1:3).^2)));
%     dat.A(n).S = blur_function(dm,M);

    % Include slice-gap
    M = model_slice_gap(M,gap{n},vs);
    
    dat.A(n).J = single(reshape(M, [1 1 1 3 3]));
end
end
%==========================================================================

%==========================================================================
% init_dat()
function dat = init_dat(Nii,mat,dm,R,D,window,gap,gapunit)
% Initialise projection matrices for super-resolution
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

if nargin < 6, window  = []; end
if nargin < 7, gap     = 0; end
if nargin < 8, gapunit = '%'; end

% Slice profile
window = get_window(window,Nii);

% Slice gap
gap = get_slice_gap(gap,Nii,gapunit);

C   = numel(Nii);
dat = struct('A',[]);
for c=1:C % Loop over channels
    
    matc = mat(:,:,min(c,size(mat,3)));
    
    if iscell(Nii(c))
        dat(c) = init_A(Nii{c},matc,dm,R{c},D{c},window{c},gap{c}); 
    else
        dat(c) = init_A(Nii(c),matc,dm,R{c},D{c},window{c},gap{c});    
    end
end    
end
%==========================================================================

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
        dmn = Nii{c}(n).dim;
        
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

%==========================================================================
% model_slice_gap()
function M = model_slice_gap(M,gap,vs)
% Modify M to model a (potential) slice-gap
% _______________________________________________________________________
%  Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

S = sqrt(sum(M(1:3,1:3).^2));
R = M(1:3,1:3)/diag(S);
S = S - gap./vs;
M = R*diag(S);
end
%==========================================================================

%==========================================================================
% noise_estimate()
function [noise,mu_val,info] = noise_estimate(Scans,K)
% Estimate average noise from a series of images
% FORMAT [noise,mu_val] = noise_estimate(Scans)
% Scans  - nifti structures or spm_vol or filenames of images
% K      - Number of Rician mixture components
% noise  - standard deviation estimate
% mu_val - expectation of more intense Rician
% info   - This struct can be used for plotting the fit as:
%              plot(info.x(:),info.p,'--',info.x(:), ...
%                   info.h/sum(info.h)/info.md,'b.', ...
%                   info.x(:),info.sp,'r');
% _________________________________________________________________________
%  Copyright (C) 2012-2019 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_noise_estimate.m 7599 2019-05-30 13:50:41Z mikael $

if isa(Scans,'nifti')
    Scans = spm_vol(Scans.fname); 
elseif iscell(Scans)
    Scans = spm_vol(Scans);
end

if nargin < 2, K = 2; end

noise  = zeros(numel(Scans),1);
mu_val = zeros(numel(Scans),1);
info   = struct('x',[],'h',[],'p',[],'sp',[],'md',[]);
for i=1:numel(Scans)
    Nii = Scans(i);
    f   = spm_read_vols(Nii);
    if spm_type(Nii.private.dat.dtype(1:(end-3)),'intt')
        f(f==max(f(:))) = 0;
        x      = 0:Nii.private.dat.scl_slope:max(f(:));
        [h,x]  = hist(f(f~=0),x);
    else
        x      = (0:1023)*(max(f(:))/1023);
        f(f==max(f(:))) = 0;
        [h,x]  = hist(f(f~=0 & isfinite(f)),x);
    end
    [mg,nu,sd,info(i)] = spm_rice_mixture_mod(double(h(:)),double(x(:)),K);
    noise(i)           = min(sd);

    if nargout>=2
        x          = -nu.^2./(2*sd.^2);
        msk        = x>-20;
        Laguerre   = exp(x(msk)/2).*((1-x(msk)).*besseli(0,-x(msk)/2)-x(msk).*besseli(1,-x(msk)/2));
        Ey         = zeros(size(sd));
        Ey( msk)   = sqrt(pi*sd(msk).^2/2).*Laguerre;
        Ey(~msk)   = nu(~msk);
        mu_val(i)  = max(Ey);
    end
end
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
end
%==========================================================================

%==========================================================================
% spm_rice_mixture_mod()
function [mg,nu,sig,info] = spm_rice_mixture_mod(h,x,K)
% Fit a mixture of Ricians to a histogram
% FORMAT [mg,nu,sig] = spm_rice_mixture(h,x,K)
% h    - histogram counts
% x    - bin positions (plot(x,h) to see the histogram)
% K    - number of Ricians
% mg   - integral under each Rician
% nu   - "mean" parameter of each Rician
% sig  - "standard deviation" parameter of each Rician
% info - This struct can be used for plotting the fit as:
%            plot(info.x(:),info.p,'--',info.x(:), ...
%                 info.h/sum(info.h)/info.md,'b.', ...
%                 info.x(:),info.sp,'r');
%
% An EM algorithm is used, which involves alternating between computing
% belonging probabilities, and then the parameters of the Ricians.
% The Koay inversion technique is used to compute the Rician parameters
% from the sample means and standard deviations. This is described at
% https://en.wikipedia.org/wiki/Rician_distribution
%__________________________________________________________________________
% Copyright (C) 2012-2019 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_rice_mixture.m 7595 2019-05-23 13:48:53Z mikael $

mg  = ones(K,1)/K;
nu  = (0:(K-1))'*max(x)/(K+1);
sig = ones(K,1)*max(x)/K/10;
lam = (sum(x.*h)/sum(h)/K).^2;

m0 = zeros(K,1);
m1 = zeros(K,1);
m2 = zeros(K,1);
ll = -Inf;
for iter=1:10000
    p  = zeros(numel(x),K);
    for k=1:K
        % Product Rule
        % p(class=k, x | mg, nu, sig) = p(class=k|mg) p(x | nu, sig, class=k)
        p(:,k) = mg(k)*ricepdf(x(:),nu(k),sig(k)^2) + eps;
    end

    % Sum Rule
    % p(x | mg, nu, sig) = \sum_k p(class=k, x | mg, nu, sig)
    sp  = sum(p,2);
    oll = ll;
    ll  = sum(log(sp).*h(:)); % Log-likelihood
    if ll-oll<1e-8*sum(h), break; end

    %fprintf('%g\n',ll);
    %md = mean(diff(x));
    %plot(x(:),p,'--',x(:),h/sum(h)/md,'b.',x(:),sp,'r'); drawnow

    % Bayes Rule
    % p(class=k | x, mg, nu, sig) = p(class=k, x | mg, nu, sig) / p(x | mg, nu, sig)
    p = bsxfun(@rdivide,p,sp);


    % Compute moments from the histograms, weighted by the responsibilities (p).
    for k=1:K
        m0(k) = sum(p(:,k).*h(:));              % Number of voxels in class k
        m1(k) = sum(p(:,k).*h(:).*x(:));        % Sum of the intensities in class k
        m2(k) = sum(p(:,k).*h(:).*x(:).*x(:));  % Sum of squares of intensities in class k
    end

    mg = m0/sum(m0); % Mixing proportions
    for k=1:K
        mu1 = m1(k)./m0(k);                                % Mean 
        mu2 = (m2(k)-m1(k)*m1(k)/m0(k)+lam*1e-3)/(m0(k)+1e-3); % Variance

        % Compute nu & sig from mean and variance
        [nu(k),sig(k)] = moments2param(mu1,mu2);
    end
    %disp([nu'; sig'])
end

if nargout >= 4
    % This info can be used for plotting the fit
    info    = struct;
    info.x  = x;    
    info.h  = h;
    info.p  = p;
    info.sp = sp;
    info.md = mean(diff(x));
end

function [nu,sig] = moments2param(mu1,mu2)
    % Rician parameter estimation (nu & sig) from mean (mu1) and variance
    % (mu2) via the Koay inversion technique.
    % This follows the scheme at
    % https://en.wikipedia.org/wiki/Rice_distribution#Parameter_estimation_.28the_Koay_inversion_technique.29
    % This Wikipedia description is based on:
    % Koay, C.G. and Basser, P. J., Analytically exact correction scheme
    % for signal extraction from noisy magnitude MR signals,
    % Journal of Magnetic Resonance, Volume 179, Issue = 2, p. 317–322, (2006)

    r     = mu1/sqrt(mu2);
    theta = sqrt(pi/(4-pi));
    if r>theta
        for i=1:256
            xi    = 2+theta^2-pi/8*exp(-theta^2/2)*((2+theta^2)*besseli(0,theta^2/4)+theta^2*besseli(1,theta^2/4))^2;
            g     = sqrt(xi*(1+r^2)-2);
            if abs(theta-g)<1e-6, break; end
            theta = g;
        end
        if ~isfinite(xi), xi = 1; end
        sig = sqrt(mu2)/sqrt(xi);
        nu  = sqrt(mu1^2+(xi-2)*sig^2);
    else
        nu  = 0;
        sig = (2^(1/2)*(mu1^2 + mu2)^(1/2))/2;
    end
end

function p = ricepdf(x,nu,sig2)
    % Rician PDF
    % p = ricepdf(x,nu,sig2)
    % https://en.wikipedia.org/wiki/Rice_distribution#Characterization
    p       = zeros(size(x));
    tmp     = -(x.^2+nu.^2)./(2*sig2);
    msk     = (tmp > -95) & (x*(nu/sig2) < 85) ; % Identify where Rice probability can be computed
    p(msk)  = (x(msk)./sig2).*exp(tmp(msk)).*besseli(0,x(msk)*(nu/sig2)); % Use Rician distribution
    p(~msk) = (1./sqrt(2*pi*sig2))*exp((-0.5/sig2)*(x(~msk)-nu).^2);      % Use Gaussian distribution
end
end
%==========================================================================