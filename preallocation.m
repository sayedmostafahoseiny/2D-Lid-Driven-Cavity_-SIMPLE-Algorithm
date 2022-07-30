function [u,uStar,uPrime,v,vStar,vPrime,p,pStar,pPrime,dU,dV,uOld,vOld]= preallocation(Nx,Ny)

% preallocate memory for u, v, p matrices
u = zeros(Nx+1,Ny+2);
uStar = zeros(Nx+1,Ny+2);
uPrime = zeros(Nx+1,Ny+2);

v = zeros(Nx+2,Ny+1);
vStar = zeros(Nx+2,Ny+1);
vPrime = zeros(Nx+2,Ny+1);

p = zeros(Nx+2,Ny+2);
pStar =zeros(Nx+2,Ny+2);
pPrime =zeros(Nx+2,Ny+2);

dU =zeros(Nx+1,Ny+2);
dV =zeros(Nx+2,Ny+1);

% initial guess
u(:,:) = 0.0;
v(:,:) = 0.0;
p(:,:) = 0.0;
uOld = u;
vOld = v;
pStar = p;

end

