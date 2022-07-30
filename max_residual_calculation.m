function [max_residual] = max_residual_calculation(Nx,Ny,u,uOld,v,vOld)

uX = u(2:Nx,2:Ny+1);
uXO = uOld(2:Nx,2:Ny+1);
vX = v(2:Nx+1,2:Ny);
vXO = vOld(2:Nx+1,2:Ny);
cmax1=max(max(abs((uX-uXO))));
cmax2=max(max(abs((vX-vXO))));
cmax = max([cmax1, cmax2]);
max_residual = cmax ;

end

