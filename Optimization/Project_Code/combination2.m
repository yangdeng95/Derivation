% Random + Newtons method

close all
clc

func = @(x) exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2))+sin(x(3)); 
options.maxIter = 10000;
options.display = 1;
bounds = [-10 10;-10 10;0 1];
tic
[x_opt,x_val] = Combine2(func,bounds,options);
time = toc;
fprintf(['The optimal value is ' num2str(x_val) '.\n'])
fprintf(['The elapsed time is ' num2str(time) '.\n'])

function [x_opt,x_val] = Combine2(func,bounds,options)
    s = size(bounds,1);
    N = 20;                % number of iterations
    X = zeros(s,N);
    X(:,1) = randomsamp(func,bounds);
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
    if options.display == 1
        fprintf(['Iteration ' num2str(N) ':the optimal value is ' num2str(x_val) '.\n']);
    end
    
% define the function of random search
    function x_r = randomsamp(func,bounds)
        x_samp = zeros(size(bounds,1),1);
        for kk = 1:size(bounds,1)
            a = bounds(kk,1);
            b = bounds(kk,2);
            x_samp(kk) = a+(b-a)*rand();
        end
        f_max = func(x_samp);
        for jj = 1:1000
            for kk = 1:size(bounds,1)
                a = bounds(kk,1);
                b = bounds(kk,2);
                x_samp(kk) = a+(b-a)*rand();
            end
            f_now = func(x_samp);
            if f_now > f_max
                f_max = f_now;
                x_mark = x_samp;
            end
        end
        x_r = x_mark;
    end
end