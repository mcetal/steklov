function [SLP_mat] = SLP_EVAL_NEW(Np,na,xs,ws,F_int,R,ds,A,B,cx,cy)
%
% Evaluates SLP with density sigma using Alpert's quadrature rules
   SLP_mat = zeros(Np,Np);
   dth = 2*pi/Np;
   numpoints = 15;
   n = Np/2;
%
% Trapezoid Rule First - exclude [i-na+1:i+na-1] in sum
   for i = 1:Np
%       disp(['i in trapezoid = ',num2str(i)])
       R_i = R(i,:);
       for k = 0:Np-2*na
          j = mod(i + na + k,Np); 
          if j==0
              j = Np;
          end
%          disp(['   j = ',num2str(j)])
          dR = R_i - R(j,:);
          Kern = log(norm(dR,2))/(2.d0*pi);
          SLP_mat(i,j) = dth*Kern*ds(j);
       end
   end 
%
% Now add in contribution from [i-na+1:i+na-1]
   for i = 1:Np
       alpha_i = (i-1)*dth;
       R_i = R(i,:);
           alpha_node = alpha_i + xs*dth;
           [xnode,ynode,dxnode,dynode,d2xnode,d2ynode,...
               dsnode] = curve_param(alpha_node,A,B,cx,cy);
       G = circshift(F_int,[0 i-1]);
       Kern = zeros(size(xs));
       for inode = 1:2*numpoints
           R_node = [xnode(inode) ynode(inode)];
           dR = R_i - R_node;
           Kern(inode) = log(norm(dR,2))/(2.d0*pi);
       end
       SLP_mat(i,:) = SLP_mat(i,:) + dth*(dsnode.*Kern.*ws)'*G;
   end

