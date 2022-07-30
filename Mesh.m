function [Nx,Ny,dx,x,dy,y] = Mesh(L1,L2)
 
% initialize grid matrix (Mesh)
Nx = 50; % number of cells in x direction
dx = L1/Nx;
x = 0:dx:L1;
Ny = 50; % number of cells in y direction
dy = L2/Ny;
y = 0:dy:L2;

end % End of function

