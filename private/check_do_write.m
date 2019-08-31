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