function [u,v,p] = var_corrections(Nx,Ny,alphaP,alphaU,p,pStar,pPrime,u,uStar,uPrime,v,vStar,vPrime)

% p corrections with under-relaxation
 for i = 2:Nx+1
 for j = 2:Ny+1
 p(i,j)= pStar(i,j) + alphaP * pPrime(i,j);
 end
 end
 % u corrections
 for i = 2:Nx
 for j = 2:Ny+1
 u(i,j)=uStar(i,j)+ alphaU * uPrime(i,j);
 end
 end
 % v corrections
 for i = 2:Nx+1
 for j = 2:Ny
 v(i,j)=vStar(i,j)+ alphaU * vPrime(i,j);
 end
 end
 
end

