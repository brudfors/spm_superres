%==========================================================================
% alloc_vars()
function [Nii_y,Nii_z,Nii_w,Nii_Dy,Nii_H] = alloc_vars(WriteTmpNii,Nii_x,dm,mat,DirOut,Verbose)
if Verbose, fprintf('Allocating niftis...'); end

if ~(exist(DirOut,'dir') == 7) 
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