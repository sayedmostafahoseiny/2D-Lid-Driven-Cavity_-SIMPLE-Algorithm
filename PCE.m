function [pPrime] = PCE(Nx,Ny,dx,dy,dU,dV,uStar,vStar,pPrime,n,maxNMiter)

 % solve pressure correction equation (PCE)
  % setup coefficients
 for i = 2:Nx+1
 for j = 2:Ny+1
 % boundary coefficients
 aPE(i,j) = dU(i,j)*dy;
 aPW(i,j) = dU(i-1,j)*dy;
 aPN(i,j) = dV(i,j)*dx;
 aPS(i,j) = dV(i,j-1)*dx;
 % central coefficient
 aPP(i,j) = aPE(i,j)+aPW(i,j)+aPN(i,j)+aPS(i,j);
 % RHS value
 bPP(i,j)= (uStar(i-1,j)-uStar(i,j))*dy+(vStar(i,j-1)-vStar(i,j))*dx;
 end
 end
 
 % fix pressure to zero at bottom left cell
 if n == 1
 i = 2;
 j = 2;
 aPP(i,j) = 1.0;
 aPE(i,j) = 0.0;
 aPW(i,j) = 0.0;
 aPN(i,j) = 0.0;
 aPS(i,j) = 0.0;
 bPP(i,j) = 0.0;
 end
 pPrime(:,:) = 0;
 for iter = maxNMiter
 for i = 2:Nx+1
 for j = 2:Ny+1
 pPrime(i,j) =(aPE(i,j)*pPrime(i+1,j)+ aPW(i,j)*pPrime(i-1,j)...
     +aPN(i,j)*pPrime(i,j+1)+aPS(i,j)*pPrime(i,j-1)+bPP(i,j))/aPP(i,j);
 end
 end
 end

end

