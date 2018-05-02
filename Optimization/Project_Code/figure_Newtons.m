% figure for Newtons method
close all
clc

func = @(x) exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2))+sin(2*pi*x(3)); 
bounds = [-10 10;-10 10;0 1];
N_list = 5:5:30;
F = zeros(1,length(N_list));


for kk = 1:3
    for jj = 1:length(N_list)
        N = N_list(jj);
        x_val = Newtons(func,bounds,N);
        F(1,jj) = x_val;
    end
    scatter(N_list,F,'g','filled');
    grid on
    hold on
end   
title('Newtons Method');
xlabel('iteration number');
ylabel('Max value');
legend('Newtons Method');

function x_val = Newtons(func,bounds,N)
    s = size(bounds,1);
    X = zeros(s,N);
    x_0 = zeros(s,1);
    for ii = 1:s
        a = bounds(ii,1);
        b = bounds(ii,2);
        x_0(ii) = a+(b-a)*rand();
    end
    X(:,1) = x_0;
% two dimensions
if s == 2
    for ii = 2:N
        syms x1 x2
        F = func([x1,x2]);
        v = [x1,x2];     
        H = hessian(F,v);
        G = gradient(F,v);
        
        x1 = X(1,ii-1);
        x2 = X(2,ii-1);
        H_value = subs(H);
        G_value = subs(G);
        M = H_value\G_value;
        X(:,ii) = X(:,ii-1)-M;
    end
    x_opt = X(:,N);
    x_val = func(x_opt);
% three dimensions
elseif s == 3
    for ii = 2:N
        syms x1 x2 x3
        F = func([x1,x2,x3]);
        v = [x1,x2,x3];     
        H = hessian(F,v);
        G = gradient(F,v);
        
        x1 = X(1,ii-1);
        x2 = X(2,ii-1);
        x3 = X(3,ii-1);
        H_value = subs(H);
        G_value = subs(G);
        M = H_value\G_value;
        X(:,ii) = X(:,ii-1)-M;
    end
    x_opt = X(:,N);
    x_val = func(x_opt);
% four dimensions
elseif s == 4
    for ii = 2:N
        syms x1 x2 x3 x4
        F = func([x1,x2,x3,x4]);
        v = [x1,x2,x3,x4];     
        H = hessian(F,v);
        G = gradient(F,v);
        
        x1 = X(1,ii-1);
        x2 = X(2,ii-1);
        x3 = X(3,ii-1);
        x4 = X(4,ii-1);
        H_value = subs(H);
        G_value = subs(G);
        M = H_value\G_value;
        X(:,ii) = X(:,ii-1)-M;
    end
    x_opt = X(:,N);
    x_val = func(x_opt);
end  
end