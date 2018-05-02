% figure for random search
close all
clc

func = @(x) exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2))+sin(2*pi*x(3)); 
bounds = [-10 10;-10 10;0 1];
N_list = 100:500:10000;
F = zeros(1,length(N_list));


for kk = 1:50
    for jj = 1:length(N_list)
        N = N_list(jj);
        x_val = randomsearch(func,bounds,N);
        F(1,jj) = x_val;
    end
    scatter(N_list,F,'r*');
    grid on
    hold on
end   
title('Random Search');
xlabel('iteration number');
ylabel('Max value');
legend('Random Search');

function x_val = randomsearch(func,bounds,N)
    x_samp = zeros(size(bounds,1),1);
    for ii = 1:size(bounds,1)
        a = bounds(ii,1);
        b = bounds(ii,2);
        x_samp(ii) = a+(b-a)*rand();
    end
    f_opt = func(x_samp);
    for jj = 1:N
        for ii = 1:size(bounds,1)
            a = bounds(ii,1);
            b = bounds(ii,2);
            x_samp(ii) = a+(b-a)*rand();
        end
        f_new = func(x_samp);
        if f_new > f_opt
            f_opt = f_new;
        end
    end
    x_val = f_opt;
end