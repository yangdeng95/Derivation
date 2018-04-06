func = @(x) tan(x);
df = @(x) 1./((cos(x))^2);
x=1.57;
h = 0:0.000001:0.0008;

for1_approx = forward1_dif(func,x,h);
for1_err = abs(df(x)-for1_approx);
loglog(h,for1_err,'r','lineWidth',2);
grid minor;
hold on;


for2_approx = forward2_dif(func,x,h);
for2_err = abs(df(x)-for2_approx);
loglog(h,for2_err,'b','lineWidth',2);
hold on;


cent1_approx = central1_dif(func,x,h);
cent1_err = abs(df(x)-cent1_approx);
loglog(h,cent1_err,'g','lineWidth',2);
hold on;


cent2_approx = central2_dif(func,x,h);
cent2_err = abs(df(x)-cent2_approx);
loglog(h,cent2_err,'y','lineWidth',2);
hold on;

comp_approx = complex_dif(func,x,h);
comp_err = abs(df(x)-comp_approx);
loglog(h,comp_err,'k','lineWidth',2);
hold on;

title('Problem 1');
xlabel('h');
ylabel('error');
legend('forward1','forward2','central1','central2','complex');

function df=forward1_dif(func,x,h)
df=(func(x+h)-func(x))./h;
end

function df = forward2_dif(func,x,h)
df = ((-11/6) * func(x) + 3 * func(x+h) - ...,
    (3/2) * func(x+2*h) + (1/3) * func(x+3*h))./h;
end

function df = central1_dif(func,x,h)
df = (func(x+h)-func(x-h))./(2*h);
end

function df = central2_dif(func,x,h)
df = ((-1/60)*func(x-3*h)+(3/20)*func(x-2*h)- ...,
    (3/4)*func(x-h)+(3/4)*func(x+h)-(3/20)*func(x+2*h) ...,
    +(1/60)*func(x+3*h))./h;
end

function df = complex_dif(func,x,h)
df =imag(func(x+i*h))./h;
end
