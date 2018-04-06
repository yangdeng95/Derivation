% Problem3 Deng Yang

func = @(x,y) exp(-0.5.*(x.^2)).*exp(-0.5.*(y.^2)).*sin(2*pi.*x).*cos(pi.*y);
fx = @(x,y) 2*pi.*exp(-x.^2/2).*exp(-y.^2/2).*cos(2*pi.*x).*cos(pi.*y) - x.*exp(-x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x);
fy = @(x,y) -pi.*exp(-x.^2/2).*exp(-y.^2/2).*sin(2*pi.*x).*sin(pi.*y)-y.*exp(-x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x);
fxx = @(x,y) x.^2.*exp(-x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x)-4*pi^2.*exp(...
             -x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x)-exp(-x.^2/2).*exp(...
             -y.^2/2).*cos(pi.*y).*sin(2*pi.*x)-4*pi.*x.*exp(-x.^2/2).*exp(-y.^2/2).*cos(2*pi.*x).*cos(pi.*y);
fyy = @(x,y) y.^2.*exp(-x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x)- pi^2.*exp(...
             -x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x)-exp(-x.^2/2).*exp(...
             -y.^2/2).*cos(pi.*y).*sin(2*pi.*x)+2*pi.*y.*exp(-x.^2/2).*exp(-y.^2/2).*sin(2*pi.*x).*sin(pi.*y);     
fxy = @(x,y) pi.*x.*exp(-x.^2/2).*exp(-y.^2/2).*sin(2*pi.*x).*sin(pi.*y)-2*pi.*y.*exp(...
             -x.^2/2).*exp(-y.^2/2).*cos(2*pi.*x).*cos(pi.*y)-2*pi^2.*exp(-x.^2/2).*exp(...
             -y.^2/2).*cos(2*pi.*x).*sin(pi.*y)+x.*y.*exp(-x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x);       
fyx = @(x,y) pi.*x.*exp(-x.^2/2).*exp(-y.^2/2).*sin(2*pi.*x).*sin(pi.*y)-2*pi.*y.*exp(...
             -x.^2/2).*exp(-y.^2/2).*cos(2*pi.*x).*cos(pi.*y)-2*pi^2.*exp(-x.^2/2).*exp(...
             -y.^2/2).*cos(2*pi.*x).*sin(pi.*y)+x.*y.*exp(-x.^2/2).*exp(-y.^2/2).*cos(pi.*y).*sin(2*pi.*x);

func_max = func(0.2375,0);
n_list = 1:1:5;
error_list = zeros(1,length(n_list));
for ii = 1:length(n_list)
    n = n_list(ii);
    y4 = newtons(func,fx,fy,fxx,fxy,fyx,fyy,0.065,0,n);
    error = abs(y4-func_max)/func_max;
    error_list(1,ii) = error;
end
% plot error-iteration
plot(n_list,error_list(1,:),'lineWidth',2);
title('Find the global maximum of the equation');
xlabel('iteration number');
ylabel('relative error');
legend('Newtons Method');

% newton's method
function sum = newtons(func,fx,fy,fxx,fxy,fyx,fyy,x,y,n)
    I = [x y]';
for ii = 1:n   
    a1 = fx(I(1),I(2));
    a2 = fy(I(1),I(2));
    b11 = fxx(I(1),I(2));
    b12 = fxy(I(1),I(2));
    b21 = fyx(I(1),I(2));
    b22 = fyy(I(1),I(2));
    H = [b11 b12;b21 b22];
    g = [a1 a2]';
    I = I - H\g;
end
    sum = func(I(1),I(2));
end
