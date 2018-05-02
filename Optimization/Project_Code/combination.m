% Random + Univarite method
close all
clc

func = @(x) exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2))+sin(x(3)); 
options.maxIter = 10000;
options.display = 1;
bounds = [-10 10;-10 10;0 1];
tic
[x_opt,x_val] = Combine1(func,bounds,options);
time = toc;
fprintf(['The optimal value is ' num2str(x_val) '.\n'])
fprintf(['The elapsed time is ' num2str(time) '.\n'])

function [x_opt,x_val] = Combine1(func,bounds,options)
% set the original point
    s = size(bounds,1);
    x_new = randomsamp(func,bounds);

    N = 2000;     % every line computational time
    h_list = zeros(1,s);
    for ii = 1:s
        h_list(ii) = (bounds(ii,2)-bounds(ii,1))/N;     % step
    end
    
    number = 0;
    f_opt = func(x_new);
    x_max = x_new;
    
    while number < 5
    for ii = 1:size(bounds,1)
        for time = 0:N
            x_new(ii) = bounds(ii,1)+h_list(ii)*time;
            f_new = func(x_new);
            if f_new > f_opt
                f_opt = f_new;
                x_max = x_new;
            end
        end
        x_new = x_max;
    end
    number = number+1;
    end
    x_val = f_opt;
    x_opt = x_max;
    if options.display == 1
        fprintf(['Iteration ' num2str(number*s+1000) ':the optimal value is ' num2str(x_val) '.\n']);
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