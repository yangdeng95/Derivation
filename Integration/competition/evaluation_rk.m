%% Challenge 3: Ball Bouncing on a Surface %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SysMat = [0 0; 1 0] ;
v0 = -5; x0 = 10 ;
ICs = [v0; x0] ; 
trange = [0, 10] ;
force = @(t,y) ballfunc(t,y); 
[true_t,true_d] = ballbounce(x0,v0,trange(2)) ;

h_list = 0.00001:0.0001:0.001;
time_list1 = zeros(1,length(h_list));
time_list2 = zeros(1,length(h_list));
time_list3 = zeros(1,length(h_list));
error_list1 = zeros(1,length(h_list));
error_list2 = zeros(1,length(h_list));
error_list3 = zeros(1,length(h_list));


 %Assess computational time and accuracy
for ii = 1:length(h_list)  
  
  h = h_list(ii);
  tic
  [t,d] = int_RK4(force,SysMat,trange,ICs,h);
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
  [t,d] = int_RK4_38(force,SysMat,trange,ICs,h);
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
  [t,d] = int_Fehlberg_5(force,SysMat,trange,ICs,h);
  time = toc;
  if size(d,1) == length(SysMat)
    d = d';
    t = t';
  end
  error = (d(:,2)-interp1(true_t',true_d',t)) ;
  error_list3(ii) = norm(error(2:end).*diff(t)) ;
  time_list3(ii) = time;
end

figure(1)
loglog(h_list,error_list1,'r','linewidth',2);
grid minor
hold on
loglog(h_list,error_list2,'b','linewidth',2);
grid minor
hold on
loglog(h_list,error_list3,'g','linewidth',2);
grid minor
hold on

title('error comparison');
xlabel('h');
ylabel('error');
legend('RungeKutta','RungeKutta3/8point','Fehlberg_5');

figure(2)
loglog(h_list,time_list1,'r','linewidth',2);
grid minor
hold on
loglog(h_list,time_list2,'b','linewidth',2);
grid minor
hold on
loglog(h_list,time_list3,'g','linewidth',2);
grid minor
hold on

title('time comparison');
xlabel('h');
ylabel('time');
legend('RungeKutta','RungeKutta3/8point','Fehlberg_5');

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
