% pure random search
close all
clc

func = @(x) exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2))+sin(x(3)); 
options.maxIter = 10000;
options.display = 1;
bounds = [-10 10;-10 10;0 1];
tic
[x_opt,x_val] = randomsearch(func,bounds,options);
time = toc;
fprintf(['The optimal value is ' num2str(x_val) '.\n'])
fprintf(['The elapsed time is ' num2str(time) '.\n'])

function [x_opt,x_val] = randomsearch(func,bounds,options)
    x_samp = zeros(size(bounds,1),1);
    for ii = 1:size(bounds,1)
        a = bounds(ii,1);
        b = bounds(ii,2);
        x_samp(ii) = a+(b-a)*rand();
    end
    f_opt = func(x_samp);
    N = options.maxIter;
    for jj = 1:N
        for ii = 1:size(bounds,1)
            a = bounds(ii,1);
            b = bounds(ii,2);
            x_samp(ii) = a+(b-a)*rand();
        end
        f_new = func(x_samp);
        if f_new > f_opt
            f_opt = f_new;
            x_mark = x_samp;
        end
    end
    x_opt = x_mark;
    x_val = f_opt;
    if options.display == 1
        fprintf(['Iteration ' num2str(jj) ':the optimal value is ' num2str(x_val) '.\n']);
    end
end