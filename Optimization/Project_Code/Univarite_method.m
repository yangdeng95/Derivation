% Univarite method
close all
clc

func = @(x) TwoD_MorletProblem(x);
options.maxIter = 10000;
options.display = 1;
bounds = [-10 10;-10 10];
tic
[x_opt,x_val] = Univarite(func,bounds,options);
time = toc;
fprintf(['The optimal value is ' num2str(x_val) '.\n'])
fprintf(['The elapsed time is ' num2str(time) '.\n'])

function [x_opt,x_val] = Univarite(func,bounds,options)
% set the original point
    s = size(bounds,1);
    x_0 = zeros(s,1);
    for ii = 1:s
        a = bounds(ii,1);
        b = bounds(ii,2);
        x_0(ii) = a+(b-a)*rand();
    end
    x_new = x_0;
    N = 2000;     % every line computational time
    h_list = zeros(1,s);
    for ii = 1:s
        h_list(ii) = (bounds(ii,2)-bounds(ii,1))/N;     % step
    end
    number = 0;
    f_opt = func(x_new);
    x_max = x_new;
% find the maximum
    while number < 3
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
        fprintf(['Iteration ' num2str(number*s) ':the optimal value is ' num2str(x_val) '.\n']);
    end
end