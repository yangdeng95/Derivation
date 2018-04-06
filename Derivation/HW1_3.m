% set the parameter
L = 1;
T = 20;
t = 0.001;

N_list = [8 16 32 64 128 256];
h_list = L./(N_list-1);
M = T./t + 1;
error_list = zeros(1,length(N_list));

% set the true value
N_sta = 256;
h_sta = L./(N_sta-1);
u_sta = diag(ones(1,N_sta-2));
u_sta = padarray(u_sta,2,'post')';
u_sta = u_sta+circshift(u_sta,1,2)*(2*(h_sta^2)/(t^2)-2) ...,
    +circshift(u_sta,2,2);
u_sta = u_sta*(t^2)/(h_sta^2);
u_sta = padarray(u_sta,1,'both');

w_sta = zeros(N_sta,M);
w_sta(N_sta,2) = sin(t);
for jj = 3:M
    w_sta(:,jj) = u_sta*w_sta(:,jj-1) - w_sta(:,jj-2);
    w_sta(N_sta,jj) = sin(t*(jj-1));
end

% compute the approximation 
for ii = 1:length(h_list)
    N = N_list(ii);
    h = h_list(ii);
    
    u = diag(ones(1,N-2));
    u = padarray(u,2,'post')';
    u = u+circshift(u,1,2)*(2*(h^2)/(t^2)-2)+circshift(u,2,2);
    u = u*(t^2)/(h^2);
    u = padarray(u,1,'both');
    
    w = zeros(N,M);
    w(N,2) = sin(t);
    for jj = 3:M
        w(:,jj) = u*w(:,jj-1) - w(:,jj-2);
        w(N,jj) = sin(t*(jj-1));
    end
    
    % extend the matrix 
    w_app = zeros(N_sta, M);
    step = (N-1)/(N_sta-1);
    for kk = 1:M
        p = 1:N;
        x = w(:,kk);
        q = 1:step:N;
        y = interp1(p,x,q,'linear');
        w_app(:,kk) = y;
    end
    
    error_norm = norm(w_app-w_sta);
    error_list(ii) = error_norm; 
end  

loglog(h_list,error_list,'b','linewidth',2);
grid minor
hold on

title('Problem 3');
xlabel('h');
ylabel('L_2 norm');
text(0.1,4,'time step=0.001s');