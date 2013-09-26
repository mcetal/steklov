function [F] = fourier_interp(xs,Np)
%
% Construct fourier interpolation matrix
   n = Np/2;
   dth = 2*pi/Np;
   F = zeros(length(xs),Np);
   for i = 1:Np
       sig = zeros(Np,1);
       sig(i) = 1;
       sig_hat = fft(sig)/Np;
       zsig = zeros(length(xs),Np);
       zsig(:,1) = sig_hat(1)*ones(size(xs));
       for kmode = 1:n
          zsig(:,kmode+1) = sig_hat(kmode+1)*exp(1i*kmode*xs*dth);
          zsig(:,Np-kmode+1) = ...
             sig_hat(Np-kmode+1)*exp(-1i*kmode*xs*dth); 
       end
       zsig(:,n+1) = zeros(size(xs));
       for j = 1:length(xs)
           hat_row = real(Np*ifft(zsig(j,:)));
           F(j,i) = hat_row(1);
       end
   end
