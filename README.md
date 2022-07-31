# 2D-Lid-Driven-Cavity_SIMPLE Algorithm
2D - incompressible - Newtonian - steady - isothermal flow in lid driven cavity using SIMPLE Algorithm
The dimensionless N-S equations are discretized and solved using finite vilume method
The flow is in horizontal plane so the gravitational acceleration is not considered
Hybrid  scheme is used for convective and diffusive fluxes. It is based on Pe number
An under-relaxation factor is used for calculation of u_new and P_new for avoiding divergence
The code has divided in several functions and all of them are called in the main code "cavity2D.m"
