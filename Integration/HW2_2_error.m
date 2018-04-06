% solve x''+100x=0, x(0)=1, x'(0)=0
% calculate the exact value of peak point
period_exact = zeros(1,47);
for ii = 1:47
    period_exact(1,ii) = ii*pi/10;
end
amplitude_exact = zeros(1,47);
for ii = 1:47
    amplitude_exact(1,ii) = (-1)^ii;
end


h_list = 0.0001:0.0001:0.01;
errorlist_period_1 = zeros(1,length(h_list));
errorlist_period_2 = zeros(1,length(h_list));
errorlist_amp_1 = zeros(1,length(h_list));
errorlist_amp_2 = zeros(1,length(h_list));

for jj = 1:length(h_list)
    q0 = [1;0];          % set initial condition
    h = h_list(jj);           % Time step
    t = 0:h:15;
    
    qstar = zeros(2,length(t));
    vstar = zeros(2,length(t));

    period_app1 = zeros(1,47);
    period_app2 = zeros(1,47);
    amplitude_app1 = zeros(1,47);
    amplitude_app2 = zeros(1,47);

% explicit scheme
    A1 = [1 h;-100*h 1];
    qstar(:,1) = q0;
    for ii = 1:length(t)-1
        qstar(:,ii+1) = A1*qstar(:,ii);
    end

% find the peak point of approximation
    kk = 1;
    for ii = 1:length(t)-2
        a1 = qstar(1,ii);
        a2 = qstar(1,ii+1);
        a3 = qstar(1,ii+2);
        if (a2-a1)*(a3-a2)<0
            amplitude_app1(1,kk) = a2;
            period_app1(1,kk) = h*(ii+1);
            kk = kk+1;
        end
     end

% implicit scheme
    v0 = [1;0];          % set initial condition
    mod = 1/(1+100*(h^2));
    A2 = [mod mod*h;-100*h*mod mod];
    vstar(:,1) = q0;
    for ii = 1:length(t)-1
        vstar(:,ii+1) = A2*vstar(:,ii);
    end

    % find the peak point of approximation
    kk = 1;
    for ii = 1:length(t)-2
        a1 = vstar(1,ii);
        a2 = vstar(1,ii+1);
        a3 = vstar(1,ii+2);
        if (a2-a1)*(a3-a2)<0
            amplitude_app2(1,kk) = a2;
            period_app2(1,kk) = h*(ii+1);
            kk = kk+1;
        end
    end

% calculate the error
    error_period_1 = norm(period_exact-period_app1);
    error_period_2 = norm(period_exact-period_app2);
    error_amp_1 = norm(amplitude_exact-amplitude_app1);
    error_amp_2 = norm(amplitude_exact-amplitude_app2);

% save the error
    errorlist_period_1(jj) = error_period_1;
    errorlist_period_2(jj) = error_period_2;
    errorlist_amp_1(jj) = error_amp_1;
    errorlist_amp_2(jj) = error_amp_2;
end

figure(1)
loglog(h_list,errorlist_period_1,'r',h_list,errorlist_period_2,'b');
grid minor
title('Period Error L_2 Norm');
xlabel('h');
ylabel('Period Error');
legend('explicit scheme','implicit scheme');

figure(2)
loglog(h_list,errorlist_amp_1,'r',h_list,errorlist_amp_2,'b');
grid minor
title('Dissipation Error L_2 Norm');
xlabel('h');
ylabel('Dissipation Error');
legend('explicit scheme','implicit scheme');

