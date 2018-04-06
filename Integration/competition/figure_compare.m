%% Challenge 3: Ball Bouncing on a Surface %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SysMat = [0 0; 1 0] ;
v0 = -5; x0 = 10 ;
ICs = [v0; x0] ; 
trange = [0, 10] ;
force = @(t,y) ballfunc(t,y); 
[true_t,true_d] = ballbounce(x0,v0,trange(2)) ;

plot(true_t,true_d);
grid minor
hold on

h_list = 0.00001:0.00001:0.0001;

 %Assess computational time and accuracy
for ii = 1:length(h_list) 
  h = h_list(ii);
  [t,d] = int_RK4_38(force,SysMat,trange,ICs,h);
  d_2 = d(2,:);
  plot(t,d_2);
  grid minor
  hold on
end
title('figure comparison');
xlabel('h');
ylabel('value');
text(6,8,'h = 0.00001:0.00001:0.0001;');

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