% sample code for function evaluation

func =  @(x) log(x);
d2f_truth = @(x) -x^(-2);
xval = 0.01;
h_list = [0.00001:0.0001:1];
error_list1 = zeros(1,length(h_list));
error_list2 = zeros(1,length(h_list));
error_list3 = zeros(1,length(h_list));
error_list4 = zeros(1,length(h_list));

for ii = 1:length(h_list)
    h = h_list(ii);
    d2f_1 = func_1(func,xval,h);
    error1 = abs(d2f_1-d2f_truth(xval))/abs(d2f_truth(xval));
    error_list1(ii) = error1;
end
for ii = 1:length(h_list)
    h = h_list(ii);
    d2f_2 = func_2(func,xval,h);
    error2 = abs(d2f_2-d2f_truth(xval))/abs(d2f_truth(xval));
    error_list2(ii) = error2;
end
for ii = 1:length(h_list)
    h = h_list(ii);
    d2f_3 = func_3(func,xval,h);
    error3 = abs(d2f_3-d2f_truth(xval))/abs(d2f_truth(xval));
    error_list3(ii) = error3;
end
for ii = 1:length(h_list)
    h = h_list(ii);
    d2f_4 = func_4(func,xval,h);
    error4 = abs(d2f_4-d2f_truth(xval))/abs(d2f_truth(xval));
    error_list4(ii) = error4;
end

% figure
loglog(h_list,error_list1,'b','lineWidth',2);
grid minor
hold on
loglog(h_list,error_list2,'y','lineWidth',2);
hold on
loglog(h_list,error_list3,'g','lineWidth',2);
hold on
loglog(h_list,error_list4,'r','lineWidth',2);
hold on

title('function=lnx x=0.01');
xlabel('h');
ylabel('error');
legend('central2','central4','central6','central8');

% central difference
function d2f = func_1(func,x,h)
    d2f = (func(x+h)-2*func(x)+func(x-h))/(h^2);
end

function d2f = func_2(func,x,h)
   d2f = ((-1/12)*func(x+2*h)+(4/3)*func(x+h)+ ...,
        (-5/2)*func(x)+(4/3)*func(x-h)+(-1/12)*func(x-2*h))/(h^2);
end

function d2f = func_3(func,x,h)
    d2f = ((1/90)*func(x+3*h)+(-3/20)*func(x+2*h)+(3/2)*func(x+h)+ ...,
   (-49/18)*func(x)+(3/2)*func(x-h)+(-3/20)*func(x-2*h) ...,
   +(1/90)*func(x-3*h))/(h^2);
end

function d2f = func_4(func,x,h)
    d2f = ((-1/560)*func(x+4*h)+(8/315)*func(x+3*h)+(-1/5)*func(x+2*h) ...,
        +(8/5)*func(x+h)+(-205/72)*func(x)+(8/5)*func(x-h) ...,
        +(-1/5)*func(x-2*h)+(8/315)*func(x-3*h)+(-1/560)*func(x-4*h))/(h^2);
end
