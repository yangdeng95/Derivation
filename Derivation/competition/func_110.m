function d2f = func_110(func,x)
% Central difference approximation
    h = 0.0001;
    d2f = ((-1/560)*func(x+4*h)+(8/315)*func(x+3*h)+(-1/5)*func(x+2*h) ...,
        +(8/5)*func(x+h)+(-205/72)*func(x)+(8/5)*func(x-h) ...,
        +(-1/5)*func(x-2*h)+(8/315)*func(x-3*h)+(-1/560)*func(x-4*h))/(h^2);
end
