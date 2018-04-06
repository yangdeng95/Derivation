%% Challenge 3: Ball Bouncing on a Surface %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SysMat = [0 0; 1 0] ;
v0 = -5; x0 = 10 ;
ICs = [v0; x0] ; 
trange = [0, 10] ;
force = @(t,y) ballfunc(t,y); 
[true_t,true_d] = ballbounce(x0,v0,trange(2)) ;

h_list = 0.00001:0.00005:0.001;
time_list1 = zeros(1,length(h_list));
time_list2 = zeros(1,length(h_list));
time_list3 = zeros(1,length(h_list));
time_list4 = zeros(1,length(h_list));
error_list1 = zeros(1,length(h_list));
error_list2 = zeros(1,length(h_list));
error_list3 = zeros(1,length(h_list));
error_list4 = zeros(1,length(h_list));

 %Assess computational time and accuracy
for ii = 1:length(h_list)  
  
  h = h_list(ii);
  tic
  [t,d] = int_euler(force,SysMat,trange,ICs,h);
  time = toc;
  if size(d,1) == length(SysMat)
    d = d';
    t = t';
  end
  error = (d(:,2)-interp1(true_t',true_d',t)) ;
  error_list1(ii) = norm(error(2:end).*diff(t)) ;
  time_list1(ii) = time;
end

for ii = 1:length(h_list)  
  
  h = h_list(ii);
  tic
  [t,d] = int_AdamsBashforth_2(force,SysMat,trange,ICs,h);
  time = toc ;
  if size(d,1) == length(SysMat)
    d = d';
    t = t';
  end
  error = (d(:,2)-interp1(true_t',true_d',t)) ;
  error_list2(ii) = norm(error(2:end).*diff(t)) ;
  time_list2(ii) = time;
end

for ii = 1:length(h_list)  
  
  h = h_list(ii);
  tic
  [t,d] = int_AdamsBashforth_4(force,SysMat,trange,ICs,h);
  time = toc;
  if size(d,1) == length(SysMat)
    d = d';
    t = t';
  end
  error = (d(:,2)-interp1(true_t',true_d',t)) ;
  error_list3(ii) = norm(error(2:end).*diff(t)) ;
  time_list3(ii) = time;
end

for ii = 1:length(h_list)  
  
  h = h_list(ii);
  tic
  [t,d] = int_RK4(force,SysMat,trange,ICs,h);
  time = toc ;
  if size(d,1) == length(SysMat)
    d = d';
    t = t';
  end
  error = (d(:,2)-interp1(true_t',true_d',t)) ;
  error_list4(ii) = norm(error(2:end).*diff(t)) ;
  time_list4(ii) = time;
end

figure(1)
loglog(h_list,error_list1,'b','linewidth',2);
grid minor
hold on
loglog(h_list,error_list2,'y','linewidth',2);
grid minor
hold on
loglog(h_list,error_list3,'g','linewidth',2);
grid minor
hold on
loglog(h_list,error_list4,'r','linewidth',2);
grid minor
hold on
title('error comparison');
xlabel('h');
ylabel('error');
legend('Euler','AdamsBashforth2','AdamsBashforth4','RungeKutta');

figure(2)
loglog(h_list,time_list1,'b','linewidth',2);
grid minor
hold on
loglog(h_list,time_list2,'y','linewidth',2);
grid minor
hold on
loglog(h_list,time_list3,'g','linewidth',2);
grid minor
hold on
loglog(h_list,time_list4,'r','linewidth',2);
grid minor
hold on
title('time comparison');
xlabel('h');
ylabel('time');
legend('Euler','AdamsBashforth2','AdamsBashforth4','RungeKutta');

function force = ballfunc(~,y)
% constants
g = 9.81 ;
k = 1e6 ;
m = 0.01 ; %i.e. a ball bearing
w = sqrt(k/m) ;
zeta = 0.05 ;

if y(2) >= 0
  force = [-g; 0];
else
  force = [-2*zeta*w*y(1)-w^2*y(2)-g;0] ;
end
end


function [t,d] = ballbounce(x0,v0,T)

% constants
g = 9.81 ;
k = 1e6 ;
m = 0.01 ; %i.e. a ball bearing
w = sqrt(k/m) ;
zeta = 0.05 ;
wd = sqrt(1-zeta^2)*w ;

t0 = 0;
t_next = findImpact(t0,x0,v0) ;
X0 = x0; V0 = v0;
t = t0;
d = X0;

while t_next < T
  % Free flight:
  ts = linspace(t0,t_next,1001) ;
  t = [t ts(2:end)]; %#ok<*AGROW>
  d = [d (-1/2*g*(ts(2:end)-t0).^2+V0*(ts(2:end)-t0)+X0)] ;
  X0 = (-1/2*g*(t_next-t0)^2+V0*(t_next-t0)+X0) ;   %#ok<NASGU>
  V0 = (-g*(t_next-t0)+V0) ;
  % Rebound:
  t0 = t_next ;
  t_free = findFree(t0,V0,m,k,zeta) ;
  ts = linspace(t0,t_free,101) ;
  t = [t ts(2:end)] ;
  A = m*g/k-sqrt(-1)*(zeta*w*(m*g/k)+V0)/wd ;
  d = [d real(A*exp((-zeta*w+sqrt(-1)*wd)*(ts(2:end)-t0))-m*g/k)] ; 
  X0 = real(A*exp((-zeta*w+sqrt(-1)*wd)*(t_free-t0))-m*g/k);  
  V0 = real(A*(-zeta*w+sqrt(-1)*wd)*exp((-zeta*w+sqrt(-1)*wd)*(t_free-t0)));
  t0 = t_free ;
  t_next = findImpact(t0,X0,V0) ;
end
% Free flight:
ts = linspace(t0,t_next,1001) ;
t = [t ts(2:end)];
d = [d (-1/2*g*(ts(2:end)-t0).^2+V0*(ts(2:end)-t0)+X0)] ;
end

function t_next = findImpact(t0,X0,V0)
g = 9.81 ;
t_next = (t0+V0/g)+sqrt(V0^2/g^2+2*X0/g) ;
end

function t_next = findFree(t0,V0,m,k,zeta)
g = 9.81 ;
w = sqrt(k/m) ;
wd = sqrt(1-zeta^2)*w ;
% A = m*g/k-sqrt(-1)*(zeta*w*(m*g/k)+V0)/wd ;
% phi = atan2(imag(A),real(A)) ;
% t_next = t0+(2*pi-2*phi)/wd ;
% t_next = t0+real(log(m*g/k/A)/(sqrt(-1)*wd-zeta*w)) ;
% if t_next < t0
%   t_next = t0+real((log(m*g/k/A)+2*pi*sqrt(-1))/(sqrt(-1)*wd-zeta*w)) ;
% end
options = optimset('Display','off','TolX',eps,'TolFun',eps);
t_next = fzero(@(t) real((g/w^2-sqrt(-1)*(g/w*zeta+V0)/sqrt(1-zeta^2)/w)*exp((-zeta*w+sqrt(-1)*w*sqrt(1-zeta^2))*(t-t0))-m*g/k), ...
  [t0+pi/wd,t0+3*pi/2/wd],options) ;
end

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
