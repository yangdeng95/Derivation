% set parameter
func = @(x) exp(x);
xval = 1;

h_list = [0:0.0001:0.001];
time_list1 = zeros(1,length(h_list));
time_list2 = zeros(1,length(h_list));
time_list3 = zeros(1,length(h_list));
time_list4 = zeros(1,length(h_list));

% time 
for ii = 1:length(h_list)
    h = h_list(ii);
% assess computational time:
tic
    for cntr = 1:1000000
        d2f = func_1(func,xval,h);
    end
time = toc/1000000;
time_list1(ii) = time;
end
for ii = 1:length(h_list)
    h = h_list(ii);
% assess computational time:
tic
    for cntr = 1:1000000
        d2f = func_2(func,xval,h);
    end
time = toc/1000000;
time_list2(ii) = time;
end
for ii = 1:length(h_list)
    h = h_list(ii);
% assess computational time:
tic
    for cntr = 1:1000000
        d2f = func_3(func,xval,h);
    end
time = toc/1000000;
time_list3(ii) = time;
end
for ii = 1:length(h_list)
    h = h_list(ii);
% assess computational time:
tic
    for cntr = 1:1000000
        d2f = func_4(func,xval,h);
    end
time = toc/1000000;
time_list4(ii) = time;
end

% figure
loglog(h_list,time_list1,'b',h_list,time_list2,'y', ...,
    h_list,time_list3,'g',h_list,time_list4,'r')
grid on;
title('function=exp(x) x=1');
xlabel('h');
ylabel('error');
legend('central2','central4','central6','central8');

% central difference
function d2f = func_1(func,x,h)
    d2f = (func(x+h)-2*func(x)+func(x-h))./(h.^2);
end
function d2f = func_2(func,x,h)
    d2f = ((-1/12)*func(x+2*h)+(4/3)*func(x+h)+ ...,
        (-5/2)*func(x)+(4/3)*func(x-h)+(-1/12)*func(x-2*h))./(h.^2);
end

function d2f = func_3(func,x,h)
    d2f = ((1/90)*func(x+3*h)+(-3/20)*func(x+2*h)+(3/2)*func(x+h)+ ...,
   (-49/18)*func(x)+(3/2)*func(x-h)+(-3/20)*func(x-2*h) ...,
   +(1/90)*func(x-3*h))./(h.^2);
end

function d2f = func_4(func,x,h)
    d2f = ((-1/560)*func(x+4*h)+(8/315)*func(x+3*h)+(-1/5)*func(x+2*h) ...,
        +(8/5)*func(x+h)+(-205/72)*func(x) +(8/5)*func(x-h) ...,
        +(-1/5)*func(x-2*h)+(8/315)*func(x-3*h)+(-1/560)*func(x-4*h))./(h.^2);
end