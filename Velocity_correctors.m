function [uPrime,vPrime] = Velocity_correctors(Nx,Ny,uPrime,dU,pPrime,vPrime,dV)

% calculate velocity corrections
% u corrections
 for i = 2:Nx
 for j = 2:Ny+1
 uPrime(i,j)=dU(i,j)*(pPrime(i,j)-pPrime(i+1,j));
 end
 end
 % v corrections
 for i = 2:Nx+1
 for j = 2:Ny
 vPrime(i,j)=dV(i,j)*(pPrime(i,j)-pPrime(i,j+1));
 end
 end
 
end

