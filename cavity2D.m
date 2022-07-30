% solving the Steady, Incompressible, Isothermal, 2D, Laminar flow in Lid driven cavity
% finite volume method + SIMPLE Algorithm with Hybrid Scheme For Convective and diffusive fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; format long

% Defining parameters
nu = input ( ' Enter the kinematic viscisity (m^2/s) = ' ) ;
[L1,L2,ulid,Re,alphaU,alphaP,maxNMiter,err_criteria,Min_Iteration,max_residual]=parameters(nu);
fprintf('\n Reynolds number = %05e \n',Re); 
disp ( ' ********************************************* ')

% Defining computational Grid (Mesh)
[Nx,Ny,dx,x,dy,y] = Mesh(L1,L2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preallocation of u,v,w,p and their star and their prime Tensors
[u,uStar,uPrime,v,vStar,vPrime,p,pStar,pPrime,dU,dV,uOld,vOld]= preallocation(Nx,Ny);

% setting Boundary conditions
[u,v] = setting_BCs(u,v,Nx,Ny,ulid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = tic; % for calculation of computational time
% SIMPLE algorithm ( SIMPLE Loop )
n = 1 ; % counter for the SIMPLE loop
while ( max_residual > err_criteria || n < Min_Iteration)
    
 % STEP 1a: solve x-momentum as uStar
[auP,dU,uStar] = x_mom(Nx,Ny,dx,dy,uOld,vOld,Re,alphaU,maxNMiter,pStar,uStar,dU);
    
 % STEP 1b: solve y-momentum as vStar
[avP,dV,vStar] = y_mom(Nx,Ny,dx,dy,uOld,vOld,Re,alphaU,maxNMiter,pStar,vStar,dV);

 % STEP 2: solve pressure correction equation (PCE)
[pPrime] = PCE(Nx,Ny,dx,dy,dU,dV,uStar,vStar,pPrime,n,maxNMiter);

 % STEP 3: calculate velocity corrections
[uPrime,vPrime] = Velocity_correctors(Nx,Ny,uPrime,dU,pPrime,vPrime,dV);

 % STEP 4: correct u, v and p
 [u,v,p] = var_corrections(Nx,Ny,alphaP,alphaU,p,pStar,pPrime,u,uStar,uPrime,v,vStar,vPrime);
 
  % check for convergence ( max residual for breaking the SIMPLE loop )
[max_residual] = max_residual_calculation(Nx,Ny,u,uOld,v,vOld);
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % STEP 5: update cell velocities ( update old values )
 uOld = u;
 vOld = v;
 pStar = p;
 n = n + 1 ; % update the counter
 
end % End of the SIMPLE loop
telapsed = toc(tstart); % for the calculation of computational time
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processing ( profiles - contours - streamlines )
postProcessing(n,telapsed,max_residual,x,y,Nx,Ny,u,v,p,dx,dy,L1,L2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









