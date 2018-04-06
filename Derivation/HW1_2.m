func = @(x) (1/6)*x^3-(1/2)*x^2;

% parameter setting
L = 1;

N_list = [10:5:200];
h_list = L./(N_list-1);

error_list = zeros(1,length(N_list));
error_tip = zeros(1,length(N_list));

for ii = 1:length(h_list)
    N = N_list(ii);
    h = h_list(ii);
    x1=diag(ones(1,N-4));
    x1=x1/(h^4);
    x=padarray(x1,4,'post')';
    x=x+circshift(x,1,2)*-4+circshift(x,2,2)*6 ...,
        +circshift(x,3,2)*-4+circshift(x,4,2);
    x=padarray(x,2,'both');
    x(1,1)=1;
    x(2,2)=1;
    x(N-1,N-2)=1/(2*(h^2));
    x(N-1,N-1)=-1/(h^2);
    x(N-1,N)=1/(2*(h^2));
    x(N,N-3)=-1/(6*(h^3));
    x(N,N-2)=1/(2*(h^3));
    x(N,N-1)=-1/(2*(h^3));
    x(N,N)=1/(6*(h^3));
    a = zeros(N,1);
    a(N,1) = 1;
    W=x\a;
    
    m = [0:h:L]';
    Y = (1/6)*m.^3-(1/2)*m.^2;
    
    error = (Y-W).^2;
    error_norm = sqrt(sum(error)*h);
    error_list(ii) = error_norm;
    
    error_tipdef = abs(func(1)-W(N,1));
    error_tip(ii) = error_tipdef;
    
end

figure(1)
loglog(h_list,error_list,'b');
grid minor
title('Problem 2 L_2 norm');
xlabel('h');
ylabel('L_2 norm');

figure(2)
loglog(h_list,error_tip,'r');
grid minor
title('Problem 2 Tip Deflection');
xlabel('h');
ylabel('Tip Deflection ErrorL_2 norm');
