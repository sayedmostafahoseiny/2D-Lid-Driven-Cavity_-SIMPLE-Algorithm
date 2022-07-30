function [avP,dV,vStar] = y_mom(Nx,Ny,dx,dy,uOld,vOld,Re,alphaU,maxNMiter,pStar,vStar,dV)

    % STEP 1b: solve y-momentum as vStar
for i = 2:Nx+1
 for j = 2:Ny
 % convective flux using central differencing
 vCe = dy*0.5*(uOld(i,j)+uOld(i,j+1));
 vCw = dy*0.5*(uOld(i-1,j)+uOld(i-1,j+1));
 vCn = dx*0.5*(vOld(i,j)+vOld(i,j+1));
 vCs = dx*0.5*(vOld(i,j)+vOld(i,j-1));
 vDx = dy/(dx*Re);
 vDy = dx/(dy*Re);
 % boundary coefficients using hybrid scheme
 avE_h = [-vCe,(vDx-0.5*vCe),0];
 avW_h = [vCw,(vDx+0.5*vCw),0];
 avN_h = [-vCn,(vDy-0.5*vCn),0];
 avS_h = [vCs,(vDy+0.5*vCs),0];
 avE(i,j) = max(avE_h);
 avW(i,j) = max(avW_h);
 avN(i,j) = max(avN_h);
 avS(i,j) = max(avS_h);
 % central coefficient using applied under-relaxation
 avP(i,j) = (avE(i,j)+avW(i,j)+avN(i,j)+avS(i,j)+(vCe-vCw)+(vCn-vCs));
 
 % set u velocity component of PCE
 dV(i,j) = dx/avP(i,j);
 end
 end
 % use previous calculation as initial guess in numerical method
 vStar = vOld;
 for iter = maxNMiter
 for i = 2:Nx+1
 for j = 2:Ny
 vStar(i,j) =(alphaU/avP(i,j))*(avE(i,j)*vStar(i+1,j)+ avW(i,j)*vStar(i-1,j)...
     + avN(i,j)*vStar(i,j+1)+avS(i,j)*vStar(i,j-1)- dx*(pStar(i,j+1)-pStar(i,j)))...
     + (1-alphaU)*vOld(i,j);
 end
 end
 end

end

