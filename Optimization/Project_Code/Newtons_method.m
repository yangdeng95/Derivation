% Newton's Method
close all
clc

func = @(x) TwoD_MorletProblem(x);
options.maxIter = 10000;
options.display = 1;
bounds = [-10 10;-10 10];
tic
[x_opt,x_val] = Newtons(func,bounds,options);
time = toc;
fprintf(['The optimal value is ' num2str(x_val) '.\n'])
fprintf(['The elapsed time is ' num2str(time) '.\n'])

function [x_opt,x_val] = Newtons(func,bounds,options)
    s = size(bounds,1);
    N = 20;                % number of iterations
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
    if options.display == 1
        fprintf(['Iteration ' num2str(N) ':the optimal value is ' num2str(x_val) '.\n']);
    end
end