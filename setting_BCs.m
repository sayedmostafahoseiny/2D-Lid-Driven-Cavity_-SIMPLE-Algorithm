function [u,v] = setting_BCs(u,v,Nx,Ny,ulid)

 % Lid Driven Cavity Boundary conditions
 % apply boundary conditions
 u(:,1) = 0.0; % bottom boundary
 v(:,1) = 0.0;
 u(:,Ny+2) = ulid; % top boundary
 v(:,Ny+1) = 0.0;
 u(1,:) = 0.0; % left boundary
 v(1,:) = 0.0;
 u(Nx+1,:) = 0.0; % right boundary
 v(Nx+2,:) = 0.0;
 
end

