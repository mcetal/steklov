function [SLP_sig] = SLP_EVAL(Np,na,xs,ws,R,ds,sigma,A,B,cx,cy)
%
% Evaluates SLP with density sigma using Alpert's quadrature rules
   SLP_sig = zeros(Np,1);
%
% Interpolate sigma to nodes
   numpoints = 15;
   n = Np/2;
   dth = 2*pi/Np;
   sig_hat = fft(sigma)/Np;
   sig_vec1 = zeros(Np,15);
   sig_vec2 = zeros(Np,15);
   for inode = 1:numpoints
      zsig(1) = sig_hat(1);
      for kmode = 1:n
          zsig(kmode+1) = sig_hat(kmode+1)*exp(1i*kmode*xs(inode)*dth);
          zsig(Np-kmode+1) = ...
             sig_hat(Np-kmode+1)*exp(-1i*kmode*xs(inode)*dth); 
      end
      zsig(n+1) = 0.d0;
      sig_vec1(:,inode) = real(Np*ifft(zsig));
      zsig(1) = sig_hat(1);
      for kmode = 1:n
          zsig(kmode+1) = sig_hat(kmode+1)*exp(-1i*kmode*xs(inode)*dth);
          zsig(Np-kmode+1) = ...
             sig_hat(Np-kmode+1)*exp(1i*kmode*xs(inode)*dth); 
      end
      zsig(n+1) = 0.d0;
      sig_vec2(:,inode) = real(Np*ifft(zsig));
   end
%
% Trapezoid Rule First - exclude [i-na+1:i+na-1] in sum
   SLP_sig = zeros(Np,1);
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
          SLP_sig(i) = SLP_sig(i) + dth*Kern*ds(j)*sigma(j);
       end
   end 
%
% Now add in contribution from [i-na+1:i+na-1]
   for i = 1:Np
       alpha_i = (i-1)*dth;
       R_i = R(i,:);
       for inode = 1:numpoints
           alpha_node = alpha_i + xs(inode)*dth;
           [xnode,ynode,dxnode,dynode,d2xnode,d2ynode,...
               dsnode] = curve_param(alpha_node,A,B,cx,cy);
           R_node = [xnode ynode];
           dR = R_i - R_node;
           Kern = log(norm(dR,2))/(2.d0*pi);
           SLP_sig(i) = SLP_sig(i) ...
               + ws(inode)*sig_vec1(i,inode)*Kern*dth*dsnode;
%
           alpha_node = alpha_i - xs(inode)*dth;
           [xnode,ynode,dxnode,dynode,d2xnode,d2ynode,...
               dsnode] = curve_param(alpha_node,A,B,cx,cy);
           R_node = [xnode ynode];
           dR = R_i - R_node;
           Kern = log(norm(dR,2))/(2.d0*pi);
           SLP_sig(i) = SLP_sig(i) ...
               + ws(inode)*sig_vec2(i,inode)*Kern*dth*dsnode;
       end
   end

