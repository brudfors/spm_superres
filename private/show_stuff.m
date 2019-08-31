%==========================================================================
% show_stuff()
function show_stuff(in,nam,np,verbose,ShowZoomed,isrgb,nr,nc,figname)
if nargin < 4, verbose    = 2; end
if nargin < 5, ShowZoomed = false; end
if nargin < 6, isrgb      = false; end
if nargin < 7, nr         = 2; end
if nargin < 8, nc         = 2; end
if nargin < 9, figname    = '(SPM) spm_superres'; end

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