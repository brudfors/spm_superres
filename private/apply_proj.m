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