function [t,X] = solvesystem_chaud496(f,g,t0,tN,x0,h)

steps = (tN-t0)/h;
x1(1) = x0(1);
x2(1) = x0(2);
t(1) = t0;

for i = 1:steps
   t(i+1) = t(i) + h;
   x1(i+1) = 0;
   x2(i+1) = 0;
   
   mx1 = f(t(i), x1(i), x2(i));
   mx2 = g(t(i), x1(i), x2(i));
   kx1 = f(t(i+1), x1(i) + mx1*h, x2(i)+ mx2*h);
   kx2 = g(t(i+1), x1(i) + mx1*h, x2(i)+ mx2*h);

   x1(i+1) = x1(i) + (h/2)*(mx1+kx1);
   x2(i+1) = x2(i) + (h/2)*(mx2+kx2);
end

X.x1 = x1;
X.x2 = x2;
end