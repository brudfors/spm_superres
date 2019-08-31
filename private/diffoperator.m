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