% figure for random-univarite

close all
clc

func = @(x) exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2))+sin(2*pi*x(3)); 
bounds = [-10 10;-10 10;0 1];
N_list = 3:3:30;
F = zeros(1,length(N_list));

for kk = 1:3
    for jj = 1:length(N_list)
        N = N_list(jj);
        x_val = Combine1(func,bounds,N);
        F(1,jj) = x_val;
    end
    scatter(N_list,F,'y','filled');
    grid on
    hold on
end   
title('Random-Univarite Method');
xlabel('iteration number');
ylabel('Max value');
legend('Random-Univarite Method');

function x_val = Combine1(func,bounds,Num)
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
    n_iter = Num/s;
    
    while number < n_iter
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