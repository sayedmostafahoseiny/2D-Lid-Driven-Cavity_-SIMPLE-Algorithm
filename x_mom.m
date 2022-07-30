function [auP,dU,uStar] = x_mom(Nx,Ny,dx,dy,uOld,vOld,Re,alphaU,maxNMiter,pStar,uStar,dU)

 % STEP 1a: solve x-momentum as uStar
 % setup coefficients
  for i = 2:Nx
 for j = 2:Ny+1
 % convective flux using central differencing
 uCe = dy*0.5*(uOld(i,j)+uOld(i+1,j));
 uCw = dy*0.5*(uOld(i,j)+uOld(i-1,j));
 uCn = dx*0.5*(vOld(i,j)+vOld(i+1,j));
 uCs = dx*0.5*(vOld(i,j-1)+vOld(i+1,j-1));
 uDx = dy/(dx*Re);
 uDy = dx/(dy*Re);
 
 % boundary coefficients using hybrid schemes
 
 auE_h = [-uCe,(uDx-0.5*uCe),0];
 auW_h = [uCw,(uDx+0.5*uCw),0];
 auN_h = [-uCn,(uDy-0.5*uCn),0];
 auS_h = [uCs,(uDy+0.5*uCs),0];
 auE(i,j) = max(auE_h);
 auW(i,j) = max(auW_h);
 auN(i,j) = max(auN_h);
 auS(i,j) = max(auS_h);
 
 % central coefficient using applied under-relaxation
 auP(i,j) = (auE(i,j)+auW(i,j)+auN(i,j)+auS(i,j)+(uCe-uCw)+(uCn-uCs));
 
 % set u velocity component of PCE
 dU(i,j) = dy/auP(i,j);
 end
 end
 % use previous calculation as initial guess in numerical method
 uStar = uOld;
 
 for iter = maxNMiter
 for i = 2:Nx
 for j = 2:Ny+1
 uStar(i,j) = (alphaU/auP(i,j)) * (auE(i,j)*uStar(i+1,j)+ auW(i,j)*uStar(i-1,j)...
 + auN(i,j)*uStar(i,j+1)+auS(i,j)*uStar(i,j-1)...
 - dy*(pStar(i+1,j)-pStar(i,j)))...
 + (1-alphaU)*uOld(i,j);
 end
 end
 end
 
end % End of function 

