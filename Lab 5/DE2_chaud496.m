function [t,y] = DE2_chaud496(t0,tN,y0,y1,h,p,q,g)
    N = ceil((tN-t0)/h);
    y(1) = y0;
    y(2)  = y0 + y1*h;
    t(1) = t0;
    t(2) = t0 + h;
    dy(2) = (y(2)-y(1))/h;
    
    for n = 2:N
        y(n+1) = 2*y(n) - y(n-1) + (-p(t(n))*dy(n) - q(t(n))*y(n) + g(t(n)))*h^2; %using equation y[n+1] = 2*y[n] - y[n-1] + y''[n]*h^2 where y''[n] is given by the DE y'' = -p(t)y' - q(t)y + g(t)
        t(n+1) = t(n) + h;
        dy(n+1) = (y(n+1) - y(n))/h; %calculating the first derivative
    end
        
end