function [x_opt,x_val] = max_110(func,bounds,options)
% Deng Yang 
% set the original point
    s = size(bounds,1);
    x_new = randomsamp(func,bounds);
% every line computational time
    N = 2000;     
    h_list = zeros(1,s);
    for ii = 1:s
        h_list(ii) = (bounds(ii,2)-bounds(ii,1))/N;     % step
    end
% set the inital value
    number = 0;
    f_opt = func(x_new);
    x_max = x_new;
    
    while number < 5
    for ii = 1:s
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