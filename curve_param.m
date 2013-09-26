function [x,y,dx,dy,d2x,d2y,ds] = curve_param(alpha,A,B,cx,cy)
      x = A*cos(alpha) + cx;
      y = B*sin(alpha) + cy;
      dx = -A*sin(alpha); 
      dy = B*cos(alpha);
      d2x = -A*cos(alpha); 
      d2y = -B*sin(alpha);
      ds = sqrt(dx.^2+dy.^2);
