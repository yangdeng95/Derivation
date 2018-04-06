% Problem 1 Deng Yang
% summmarize your results in a table showing number of iterations relative
% plot the percent error for each method as a function of iteration number

% set inital number and funciton 
func = @(x) x^8-1;
funcp = @(x) 8*x^7;
n_list = 1:1:10;
x_exact = 1;
error_list = zeros(6,length(n_list));

for kk = 1:length(n_list)
    x_1 = Bisection(func,0,1.5,kk);
    error = abs(x_1-x_exact)/x_exact;
    error_list(1,kk) = error;
    x_2 = LinearInterp(func,0,1.5,kk);
    error = abs(x_2-x_exact)/x_exact;
    error_list(2,kk) = error;
    x_3 = NewtonRaphson(func,funcp,0.7,kk);
    error = abs(x_3-x_exact)/x_exact;
    error_list(3,kk) = error;
    x_4 = NewtonRaphson(func,funcp,1.5,kk);
    error = abs(x_4-x_exact)/x_exact;
    error_list(4,kk) = error;
    x_5 = Secant(func,0,1.1,kk);
    error = abs(x_5-x_exact)/x_exact;
    error_list(5,kk) = error;
    x_6 = Secant(func,1.5,0.75,kk);
    error = abs(x_6-x_exact)/x_exact;
    error_list(6,kk) = error;
end
% plot error-iteration
for jj = 1:6
    plot(n_list,error_list(jj,:),'lineWidth',2);
    hold on
end
title('Find the zero crossing');
xlabel('iteration number');
ylabel('relative error');
legend('Bisection Method','Linear Interpolation Method', ...,
    'Newton Raphson Method-lower','Newton Raphson Method-upper', ...,
    'Secant Method-lower','Secant Method-upper');

% Bisection Method
function x = Bisection(func,x_l,x_u,n)
for ii = 1:n
    x_j = 0.5*(x_l+x_u);
    yl = func(x_l);
    y2 = func(x_j);
    a = yl*y2;
    
% give the x_j to x_l or x_u
    if a < 0
        x_u = x_j;
    elseif a > 0
        x_l = x_j;
    elseif a == 0
        x = x_j;
        fprintf('%d',ii);
        return 
    end 
end   
    x = x_j;
end

% Linear Interpolation Method
function x = LinearInterp(func,x_l,x_u,n)
    y_l = func(x_l);
    y_u = func(x_u);
for ii = 1:n    
    x_j = x_l - y_l*(x_u-x_l)/(y_u-y_l);
    y_j = func(x_j);
    a = y_l*y_j;
    
% give the x_j to x_l or x_u
    if a < 0
       x_u = x_j;
    elseif a > 0
       x_l = x_j;
    elseif a == 0
       x = x_j;
       fprintf('%d',ii);
       return
    end    
end   
    x = x_j;
end

% Newton-Raphson Method
function x = NewtonRaphson(func,funcp,x_i,n)
    x_j = x_i;
for ii = 1:n
    x_j = x_j-func(x_j)/funcp(x_j);
end
    x = x_j;
end

% Secant Method
function x = Secant(func,x_1,x_2,n)
    x(1) = x_1;
    x(2) = x_2;
for ii = 3:n+2
    x(ii) = x(ii-1)-(func(x(ii-1))*(x(ii-1)-x(ii-2)))/(func(x(ii-1))-func(x(ii-2)));
end
    x = x(n+2);
end