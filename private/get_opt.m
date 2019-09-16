%==========================================================================
function opt = get_opt(opt)
% scaling of regularisation parameter (lambda)
if ~isfield(opt,'LamScl'),      opt.LamScl      = 10;           end  
% scaling of step-size parameter (rho)
if ~isfield(opt,'RhoScl'),      opt.RhoScl      = 1;            end  
% Max number of iterations
if ~isfield(opt,'MaxNiter'),    opt.MaxNiter    = 50;          end  
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
% Different test-cases: 
%   0. No testing
%   1. brainweb (superres)
%   2. brainweb (den)
%   3. lena
%   4. qmri
if ~isfield(opt,'TestCase'),    opt.TestCase    = 1;            end  
% Memory limit for when to allocate variables as niftis
if ~isfield(opt,'MaxMem'),      opt.MaxMem      = 2048;         end 
% Reconstruction voxel size, if 0, set to smallest available
if ~isfield(opt,'VoxSize'),     opt.VoxSize     = 1;            end 
% Downsample inplane resolution to 1 mm
if ~isfield(opt,'Inplane1mm'),  opt.Inplane1mm  = true;         end 
% Do just denoising, without super-resolving
if ~isfield(opt,'Denoise'),     opt.Denoise     = false;        end 
end
%==========================================================================