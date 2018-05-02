close all
clear all
clc

func = @(x) TwoD_MorletProblem(x);
options.display = 1;
bounds = [-10 10;-10 10];
tic
[x_opt,x_val] = max_020(func,bounds,options);
time = toc;
fprintf(['The optimal value is ' num2str(x_val) '.\n'])
fprintf(['The elapsed time is ' num2str(time) '.\n'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_opt,x_val] = max_020(func,bounds,options)
x_0 = zeros(size(bounds,1),1);
for cntr = 1:size(bounds,1)
    x_0(cntr) = mean(bounds(cntr,:));
end
if isfield(options,'maxIter')
    if options.display == 'iter'
        opts.Display = 'iter';
    else
        opts.Display = 'off';
    end
else
    opts.Display = 'off';
end
opts.TolX = eps;
opts.TolFun = eps;
fminfunc = @(x) 0-func(x);
[x_opt,x_val] = fminsearch(fminfunc,x_0,opts);
x_val = 0-x_val;
end