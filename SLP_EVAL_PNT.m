function [SLP_val]=SLP_EVAL_PNT(Np,R,ds,sigma,R_point)
%
% Trapezoid Rule
   SLP_val = 0;
   dth = 2*pi/Np;
   for i = 1:Np
       dR = R_point - R(i,:);
       Kern = log(norm(dR,2))/(2.d0*pi);
       SLP_val = SLP_val + dth*Kern*ds(i)*sigma(i);
   end 
