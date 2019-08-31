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