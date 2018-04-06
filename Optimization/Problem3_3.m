% Problem 3_3 Deng Yang

format long e
func = @(x,y) exp(-0.5.*(x.^2)).*exp(-0.5.*(y.^2)).*sin(2*pi.*x).*cos(pi.*y);
dx = 0.01;
x = -5:dx:5;
y = -5:dx:5;
m = 0;
n = 0;
f = func(m,n);
f_now = f;
f_exact = func(0.24,0);
 
num = 1;
error_list(num) = abs(f(num)-f_exact);
while error_list(num) > 1e-3
    ii = 1;
    f(ii) = func(m,n) ;
    while m <= 5
        error_list(num) = abs(f(ii)-f_exact)/f_exact ;
        if f(ii) > f_now
            f_now = f(ii) ;
        end
        fprintf('n=%d   f_now=%d   error=%d\n',num,f_now,error_list(num))
        if f_now == f_exact
            break
        end
        num = num+1 ;
        m = m+dx ;
        ii = ii+1 ;
        f(ii) = func(m,n) ;
    end
call = 1 ;
f(ii) = func(m,n) ;
fn(call) = func(m,n) ;
m = x(find(x == 0)+find(f == max(f))-1) ;
n = n+dx ;
 
    while n <=5
        error_list(num) = abs(f(ii)-f_exact)/f_exact ;
        if f(ii) > f_now
            f_now = f(ii) ;
        end
        fprintf('n=%d   error=%d\n',num,error_list(num))
        if f_now == f_exact
            break
        end
        num = num+1 ;
        n = n+dx ;
        ii = ii+1 ;
        f(ii) = func(m,n) ;
        call = call+1 ;
        fn(call) = func(m,n) ;
    end
        n = y(find(y == 0)+find(fn == max(fn))-1) ;
        num = num+1 ;
        ii = ii+1 ;
        n = n(:,1) ;
        f(ii) = func(m,n) ;
        error_list(num) = abs(f(ii)-f_exact) ;
        f = [] ;
        if f_now == f_exact
            break
        end
end
plot (1:num, error_list,'lineWidth',2);
title('Find the global maximum of the equation');
xlabel('iteration number');
ylabel('relative error');
legend('Univariate Method');
