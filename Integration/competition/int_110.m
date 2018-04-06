function [time,displacement] = int_110(force,SysMat,trange,ICs)
    h = 0.0001;
    t0 = trange(1);
    t1 = trange(2);
    time = t0:h:t1;
    
    displacement = zeros(length(ICs),length(time));
    displacement(:,1) = ICs;
    f = @(t,y) SysMat*y+force(t,y);
    for ii = 1:length(time)-1
        tstar = t0+(ii-1)*h;
        ystar = displacement(:,ii);
        
        k1 = f(tstar,ystar);
        k2 = f(tstar+h/3,ystar+h*(k1/3));
        k3 = f(tstar+h*2/3,ystar+h*(-1/3*k1+k2));
        k4 = f(tstar+h,ystar+h*(k1-k2+k3));
        displacement(:,ii+1) = displacement(:,ii)+h/8*(k1+3*k2+3*k3+k4);
    end
end