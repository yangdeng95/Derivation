% calculate the true value
func = @(x) (exp(x)-exp(-x))/2;
True_value = exp(10)/2+exp(-10)/2-1;
a = 0; b = 10;
N_list = 10000:-10:10;
h_list = (b-a)./N_list;

error_list1 = zeros(1,length(N_list));
error_list2 = zeros(1,length(N_list));
error_list3 = zeros(1,length(N_list));
error_list4 = zeros(1,length(N_list));
error_list5 = zeros(1,length(N_list));

% calculate error  print figure
for kk = 1:length(N_list)
    n = N_list(kk);
    Ap = Trapezoidal(func,a,b,n);
    error = abs(Ap - True_value);
    error_list1(kk) = error;
end

for kk = 1:length(N_list)
    n = N_list(kk);
    Ap = Simpson_Mid(func,a,b,n);
    error = abs(Ap - True_value);
    error_list2(kk) = error;
end

for kk = 1:length(N_list)
    n = N_list(kk);
    Ap = Simpson_38(func,a,b,n);
    error = abs(Ap - True_value);
    error_list3(kk) = error;
end

for kk = 1:length(N_list)
    n = N_list(kk);
    Ap = Gauss_2(func,a,b,n);
    error = abs(Ap - True_value);
    error_list4(kk) = error;
end

for kk = 1:length(N_list)
    n = N_list(kk);
    Ap = Gauss_legendre_4(func,a,b,n);
    error = abs(Ap - True_value);
    error_list5(kk) = error;
end

% print figure
loglog(h_list,error_list1,'b','lineWidth',2);
hold on
loglog(h_list,error_list2,'g','lineWidth',2);
hold on
loglog(h_list,error_list3,'y','lineWidth',2);
hold on
loglog(h_list,error_list4,'m','lineWidth',2);
hold on
loglog(h_list,error_list5,'r','lineWidth',2);
grid minor
hold on

title('Different Method of Integration');
xlabel('h');
ylabel('error');
legend('Trapezoidal','Simpson Midpoint', ...,
    'Simpson 3/8','Gauss 2point','Gauss legendre 4point');

% set functions
function If = Trapezoidal(func,a,b,n)
    h = (b-a)/n;
    sum = 0;
    for ii = 1:n
        sum = sum + func(a+h*(ii-1))+func(a+h*ii);
    end
    If = h*sum/2;
end

function If = Simpson_Mid(func,a,b,n)
    h = (b-a)/n;
    sum = 0;
    for ii = 1:n
        sum = sum+func(a+(ii-1)*h)+4*func(a+(ii-0.5)*h)+func(a+ii*h);
    end
    If = h*sum/6;
end

function If = Simpson_38(func,a,b,n)
    h = (b-a)/n;
    sum = 0;
    for ii = 1:n
       sum = sum + func(a+h*(ii-1))+3*func(a+h*(ii-2/3)) ...,
           +3*func(a+h*(ii-1/3))+func(a+h*ii);
    end
    If = h*sum/8;
end

function If = Gauss_2(func,a,b,n)
    h = (b-a)/n;
    If = 0;
    for ii = 1:n
       first = a+h*(ii-1);
       last = a+h*ii;
       x_1 = (first+last)/2 - h/(2*sqrt(3));
       x_2 = (first+last)/2 + h/(2*sqrt(3));
       If = If +(h/2)*(func(x_1)+func(x_2));
    end
end

function If = Gauss_legendre_4(func,a,b,n)
    c1 = 0.3478548; c2 = 0.6521452;
    x1 = -0.861136312; x2 = -0.339981044;
    x3 = 0.339981044; x4 = 0.861136312;
    h = (b-a)/n;
    If = 0; sum = 0;
    for ii = 1:n
       first = a+h*(ii-1);
       last = a+h*ii;
       R = (first+last)/2;
       sum = c1*func(h/2*x1+R)+c2*func(h/2*x2+R)+c2*func(h/2*x3+R)+c1*func(h/2*x4+R);
       If = If+sum*h/2;
    end
end