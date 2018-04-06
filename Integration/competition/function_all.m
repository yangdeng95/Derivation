function [time,displacement] = int_euler(force,SysMat,trange,ICs,h)
    t0 = trange(1);
    t1 = trange(2);
    time = t0:h:t1;
    
    displacement = zeros(length(ICs),length(time));
    displacement(:,1) = ICs;
    f = @(t,y) SysMat*y+force(t,y);
    
    for ii = 1:length(time)-1
        tstar = t0+(ii-1)*h;
        ystar = displacement(:,ii);
        
       	k = f(tstar,ystar);
        displacement(:,ii+1) = displacement(:,ii)+h*k;
    end
end

function [time,displacement] = int_AdamsBashforth_2(force,SysMat,trange,ICs,h)
    t0 = trange(1);
    t1 = trange(2);
    time = t0:h:t1;
    
    displacement = zeros(length(ICs),length(time));
    displacement(:,1) = ICs;
    f = @(t,y) SysMat*y+force(t,y);
    displacement(:,2) = displacement(:,1)+h*f(t0,displacement(:,1));
    
    for ii = 1:length(time)-2
        t_m = t0+(ii-1)*h;
        y_m = displacement(:,ii);
        t_n = t0+ii*h;
        y_n = displacement(:,ii+1);
        displacement(:,ii+2) = displacement(:,ii+1)+h*(1.5*f(t_n,y_n)-0.5*f(t_m,y_m));
    end
end

function [time,displacement] = int_AdamsBashforth_4(force,SysMat,trange,ICs,h)
    t0 = trange(1);
    t1 = trange(2);
    time = t0:h:t1;
    
    displacement = zeros(length(ICs),length(time));
    displacement(:,1) = ICs;
    f = @(t,y) SysMat*y+force(t,y);
    displacement(:,2) = displacement(:,1)+h*f(t0,displacement(:,1));
    displacement(:,3) = displacement(:,2)+ ...,
        h*(1.5*f(t0+h,displacement(:,2))-0.5*f(t0,displacement(:,1)));
    displacement(:,4) = displacement(:,3)+ ...
        h*(23/12*f(t0+2*h,displacement(:,3))-4/3*f(t0+h,displacement(:,2)) ...,
        +5/12*f(t0,displacement(:,1)));
    for ii = 1:length(time)-4
        t_1 = t0+(ii-1)*h;
        y_1 = displacement(:,ii);
        t_2 = t0+ii*h;
        y_2 = displacement(:,ii+1);
        t_3 = t0+(ii+1)*h;
        y_3 = displacement(:,ii+2);
        t_4 = t0+(ii+2)*h;
        y_4 = displacement(:,ii+3);
        displacement(:,ii+4) = y_4+h*(55/24*f(t_4,y_4)-59/24*f(t_3,y_3)+37/24*f(t_2,y_2)-3/8*f(t_1,y_1));
    end
end

function [time,displacement] = int_RK4(force,SysMat,trange,ICs,h)
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
        k2 = f(tstar+h/2,ystar+h*k1/2);
        k3 = f(tstar+h/2,ystar+h*k2/2);
        k4 = f(tstar+h,ystar+h*k3);
        displacement(:,ii+1) = displacement(:,ii)+h/6*(k1+2*k2+2*k3+k4);
    end
end

function [time,displacement] = int_RK4_38(force,SysMat,trange,ICs,h)
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

function [time,displacement] = int_Fehlberg_5(force,SysMat,trange,ICs,h)
    t0 = trange(1);
    t1 = trange(2);
    time = t0:h:t1;
   
    displacement = zeros(length(ICs),length(time));
    displacement(:,1) = ICs;
    f = @(t,y) SysMat*y+force(t,y);
    for ii = 1:length(time)-1
        ta = t0+(ii-1)*h;
        ya = displacement(:,ii);
% Calculating k1, k2, k3, k4, K5 and K6
        k1 = f(ta,ya) ;
        k2 = f(ta+(1/4)*h,ya+(1/4)*k1*h) ;
        k3 = f(ta+(3/8)*h,ya+(3/32)*h*(k1+3*k2)) ;
        k4 = f(ta+(12/13)*h,ya+(12/2197)*h*(161*k1-600*k2+608*k3)) ;
        k5 = f(ta+h,ya+(1/4104)*h*(8341*k1-32832*k2+29440*k3-845*k4)) ;
        k6 = f(ta+(0.5)*h,ya+h*(-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5)) ;
% Using 6th Order Runge-Kutta formula
    displacement(:,ii+1) = displacement(:,ii)+1/5*((16/27)*k1+(6656/2565)*k3+(28561/11286)*k4-(9/10)*k5+(2/11)*k6)*h;
   end
end