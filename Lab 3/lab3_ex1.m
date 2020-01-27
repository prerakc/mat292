function [yn,tn] = lab3_ex1(f,t0,tN,y0,h)
    tn(1) = t0;
    yn(1) = y0;
    num_steps= (tN-t0)/h;
    for i = 1:num_steps;
        yn(i+1)  = yn(i) + ((f(tn(i),yn(i)) + f(tn(i)+h, yn(i)+h*f(tn(i),yn(i))))  *(h/2));
        tn(i+1) = tn(i) + h;
    end
end