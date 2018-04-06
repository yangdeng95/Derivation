% Problem3 Deng Yang
% find the global maximum of the equation
func= @(x,y) exp(-0.5*(x^2))*exp(-0.5*(y^2))*sin(2*pi*x)*cos(pi*y);

func_max = func(0.2375,0);
n_list = 10:50:500;
error_list = zeros(2,length(n_list));
for kk = 1:50
    for ii = 1:length(n_list)
        n = n_list(ii);
        y1 = randomsearch(func,[-5,5],[-5,5],n);
        error = abs(y1-func_max)/func_max;
        error_list(1,ii) = error;
        y2 = Latinhyper(func,[-5,-5],[5,5],n);
        error = abs(y2-func_max)/func_max;
        error_list(2,ii) = error;
    end
        scatter(n_list,error_list(1,:),'r*');
        hold on
        scatter(n_list,error_list(2,:),'k');
        hold on
end
title('Find the global maximum of the equation');
xlabel('iteration number');
ylabel('relative error');
legend('Random Search','Latin Hyber Cube Sampling');
   

% random search
function sum = randomsearch(func,x_range,y_range,n)
    x1 = x_range(1);
    x2 = x_range(2);
    y1 = y_range(1);
    y2 = y_range(2);
    f_opt = func(x1,x2);
    for ii = 1:n
        x = x1+(x2-x1)*rand();
        y = y1+(y2-y1)*rand();
        f_new = func(x,y);
        if f_new > f_opt
            f_opt = f_new;
        end
    end
    sum = f_opt;
end

% latin hyper cube sampling
function f_optimal = Latinhyper(func,min,max,n)
    range=max-min;
    offset=min;
    Range=ones(n,2);
    Offset=ones(n,2);
for ii = 1:2
    Range(:,ii)=ones(n,1).*range(ii);
    Offset(:,ii)=ones(n,1).*offset(ii);
end
% returns an n-by-2 matrix, X, containing a latin hypercube sample of n
% within 0/1
    S = zeros(n,2);
for ii = 1:2
  S(:,ii) = ((randperm(n)-1+rand(1,n)))'/n;
end  
% returns an n-by-2 matrix, X, containing a latin hypercube sample of n
% corresponding to range
    LHC = Range.*S+Offset;
    f_max = func(LHC(1,1),LHC(1,2));
for jj = 2:n
    f_new = func(LHC(jj,1),LHC(jj,2));
    if f_new > f_max
        f_max = f_new;
    end
end
    f_optimal = f_max;
end
