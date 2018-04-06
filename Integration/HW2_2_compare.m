% solve x''+100x=0, x(0)=1, x'(0)=0

q0 = [1;0];        % set initial condition
h = 0.001;           % Time step
t = 0:h:15;

x_exact = cos(10*t);
qstar = zeros(2,length(t));
vstar = zeros(2,length(t));

% explicit scheme
A1 = [1 h;-100*h 1];
qstar(:,1) = q0;
for i = 1:length(t)-1
    qstar(:,i+1) = A1*qstar(:,i);
end

% implicit scheme
v0 = [1;0];          % set initial condition
mod = 1/(1+100*(h^2));
A2 = [mod mod*h;-100*h*mod mod];
vstar(:,1) = q0;
for i = 1:length(t)-1
    vstar(:,i+1) = A2*vstar(:,i);
end

plot(t,x_exact,t,qstar(1,:),t,vstar(1,:));
legend('Exact','explicit','implicit');
hold off