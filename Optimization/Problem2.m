% Problem 2 Deng Yang
% Find the maximum of the equation

% calculate the derivative
func = @(t) exp(-0.5*(t^3))*sin(pi*t);
funcp = @(t) exp(-0.5*(t^3))*(-1.5*(t^2)*sin(pi*t)+pi*cos(pi*t));
funcpp = @(t) exp(-0.5*(t^3))*((9/4*(t^4)-3*t-pi^2)*sin(pi*t)-3*(t^2)*pi*cos(pi*t));

% set the initial number
n_list = 1:1:5;
y_exact = func(0.467);
error_list = zeros(4,length(n_list));
% calculate error
for kk = 1:length(n_list)
    y1 = goldsection(func,0,1,kk);
    error = abs(y1-y_exact)/y_exact;
    error_list(1,kk) = error;
    y2 = quadraticInterp(func,0,1,0.2,kk);
    error = abs(y2-y_exact)/y_exact;
    error_list(2,kk) = error;
    y3 = Newtons(func,funcp,funcpp,0.15,kk);    % 0 and 1 both lead to the second derivation to 0
    error = abs(y3-y_exact)/y_exact;
    error_list(3,kk) = error;
    y4 = Newtons(func,funcp,funcpp,0.72,kk);
    error = abs(y4-y_exact)/y_exact;
    error_list(4,kk) = error;
end
% plot error-iteration
for jj = 1:4
    plot(n_list,error_list(jj,:),'lineWidth',2);
    hold on
end
title('Find the maximum of the equation');
xlabel('iteration number');
ylabel('relative error');
legend('Golden Section Search','Quadratic Interpolation', ...,
    'Newtons Method-lower','Newtons Method-upper');

% Golden Section Search
function y =  goldsection(func,x_l,x_u,n)
   d = (sqrt(5)-1)/2*(x_u-x_l);
   x_1 = x_l+d;
   x_2 = x_u-d;
   x_max = x_1;
for ii = 1:n
    if func(x_1) > func(x_2)
        x_max = x_1;
        x_l = x_2;
        d = (sqrt(5)-1)/2*(x_u-x_l);
        x_2 = x_1;
        x_1 = x_l+d;   
    elseif func(x_1) < func(x_2)
        x_max = x_2;
        x_u = x_1;
        d = (sqrt(5)-1)/2*(x_u-x_l);
        x_1 = x_2;
        x_2 = x_u-d;
     end
end
   y = func(x_max);
end

% Quadratic Interpolation
function y = quadraticInterp(func,x_0,x_2,x_1,n)
for ii = 1:n
    f0 = func(x_0);
    f1 = func(x_1);
    f2 = func(x_2);
    x_3 = (func(x_0)*(x_1^2-x_2^2)+func(x_1)*(x_2^2-x_0^2)+func(x_2)*(x_0^2-x_1^2))/ ...,
        (2*func(x_0)*(x_1-x_2)+2*func(x_1)*(x_2-x_0)+2*func(x_2)*(x_0-x_1));
    f3 = func(x_3);
    if f3 > f1    
    x_0 = x_1;     
    x_1 = x_3;    
    elseif f1 > f3 && f3 > f2 && x_3 < x_2 && x_3 > x_1 
            x_2 = x_3;    
    elseif f1 > f3 && f3 > f0 && x_3 > x_0 && x_3 < x_1  
            x_0 = x_3;  
    end
end  
    y = func(x_3);
end

% Newton Method
function y = Newtons(func,funcp,funcpp,x,n)
    x(1) = x;
for ii = 2:n
    x(ii) = x(ii-1)-funcp(x(ii-1))/funcpp(x(ii-1));
end
    y = func(x(n));
end
