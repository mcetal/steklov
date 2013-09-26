%
% Solves the Steklov eigenvalue problem:
%
%     lap(u) = 0 in D
%     u_n + \lambda u = 0 on \Gamma
%
% Represent u using SLP
% For now, domain is inside an ellipse
%
   close all; 
   clear;
%
% Define major and minor axis of ellipse, and centre
   A = 2; B = 1/A; cx = 0; cy = 0;
%
% Define closed curve, parametrized by alpha
   Np = 256; dth = 2*pi/Np; alpha = 0:dth:2*pi-dth;
   [x,y,dx,dy,d2x,d2y,ds] = curve_param(alpha,A,B,cx,cy);
%
% Plot
   figure(1)
   plot(x,y,'k')
   hold on
   axis equal
%
% Calculate normal and tangent vectors
   R = [x;y]';
   dR_dalpha = [dx;dy]';
   d2R_dalpha2 = [d2x;d2y]';
   N = [dy./ds;-dx./ds]';
   quiver(R(:,1),R(:,2),dR_dalpha(:,1),dR_dalpha(:,2),0.75)
   quiver(R(:,1),R(:,2),N(:,1),N(:,2),0.75)
%
% Construct discrete integral operator K
   KN = zeros(Np,Np); 
   for i = 1:Np
       R_i = R(i,:);
       N_i = N(i,:);
       kappa = norm(...
               cross([dR_dalpha(i,:) 0],[d2R_dalpha2(i,:) 0])/(ds(i))^3,2);
       for j = 1:Np
          dR = R_i - R(j,:);
          N_j = N(j,:);
          gradU = dR/(norm(dR,2)^2);
%
% Kernel for Neumann
% Modify kernel to remove rank deficiency
          KN(i,j) = dth*(dot(gradU,N_i)+1)*ds(j)/(2.d0*pi);
       end
       KN(i,i) = -0.5+dth*kappa*ds(i)/(4*pi) + dth*ds(i)/(2.d0*pi); 
   end 
%
   figure(4)
   mesh(KN)
   disp(['Condition Number = ',num2str(cond(KN))])
%
% Calculate SLP with density sigma using Alpert's quadrature rule of 
% order 16
   [xs,ws,na]=alpertquad()
   xs_all = [-flipud(xs);xs];
   ws_all = [flipud(ws);ws];
%
% Construct Fourier interpolation matrix that interpolates a vector
% to the nodes
   F_int = fourier_interp(xs_all,Np);
% Construct matrix for SLP evaluation
   SLP_mat = SLP_EVAL_NEW(Np,na,xs_all,ws_all,F_int,R,ds,A,B,cx,cy);
%
% Calculate Steklov eigenvalues
   [V,D] = eig(KN\SLP_mat);
   D = diag(D);
   disp(['Max imaginary part of eigenvalues = ',num2str(max(abs(imag(D))))])
   [rho,I] = sort(real(D),'descend');
%
% Pick a max number of eigenvalues to look at
   N_lambda = 10;
   lambda = -1./rho(1:N_lambda)
% 
% Plot solution
   NX = 100; NY = 100; dx = 2*A/NX; dy = 2*B/NY;
   [xgrid,ygrid]=meshgrid(-A:dx:A,-B:dy:B);
   eps = 0.01;
   ugrid = zeros(size(xgrid));
   for i_lambda = 1:N_lambda
       figure(1)
       clf
       sigma = V(:,I(i_lambda));
       for i =1:NX+1
          for j = 1:NY+1
             R_point = [xgrid(i,j) ygrid(i,j)];
             check = (R_point(1)/A)^2 + (R_point(2)/B)^2;
             if check<1-eps
                 ugrid(i,j) = SLP_EVAL_PNT(Np,R,ds,sigma,R_point);
             else
                 ugrid(i,j) = 0.;
             end
          end
      end
      contour(xgrid,ygrid,ugrid,30)
      axis equal
      t = ['\lambda = ',num2str(lambda(i_lambda))];
      title(t)
      pause
   end
           