% calculate y1 by y0
A = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 2 0;
     0 0 0 0 0 4];
B = [0 0 0 -1 0 0;
      0 0 0 0 -1 0;
      0 0 0 0 0 -1;
      10 -10 0 0 0 0;
      -10 30 -20 0 0 0;
      0 -20 50 0 0 0];
C = [0 0 0 1 0 0]';
  
h_list = 0.001:0.01:1;
t_max = 0:0.001:15;
M = zeros(length(h_list),length(t_max));
error_list1 = zeros(1,length(h_list));
error_list2 = zeros(1,length(h_list));
error_list3 = zeros(1,length(h_list));

% trapezoidal rule
for jj = 1:length(h_list)
   h = h_list(jj);
   t_list = 0:h:15;
   Y1 = zeros(6,length(t_list)); 
   Y1(:,1) = [1 0.5 0.25 0 0 -1]';
 
   for ii = 2:length(t_list)
      t1 = t_list(ii);
      t2 = t_list(ii-1);
      Y1(:,ii) = inv(eye(6)-0.5*h*inv(A)*(-B))*(Y1(:,ii-1)+ ...,
           0.5*h*inv(A)*(-B)*Y1(:,ii-1)+0.5*h*inv(A)*(sin(2*pi*t1) ...,
           +sin(2*pi*t2))*C);
    end 
% extend the matrix
    
   Y = interp1(t_list,Y1(1,:),t_max,'linear','extrap');
   M(jj,:) = Y;
    
end
% calculate the error L2 norm
for kk = 2:length(h_list)
   error_list1(1,kk-1) = norm(M(kk,:)-M(1,:));
end

% midpoint method
for jj = 1:length(h_list)
   h = h_list(jj);
   t_list = 0:h:15;
   Y2 = zeros(6,length(t_list));
   Y2(:,1) = [1 0.5 0.25 0 0 -1]';
   for ii = 2:length(t_list)
      t1 = t_list(ii);
      Y2(:,ii) = inv(eye(6)-0.5*h*inv(A)*(-B))*(Y2(:,ii-1)+ ...,
           0.5*h*inv(A)*(-B)*Y2(:,ii-1)+h*inv(A)*(sin(2*pi*(t1-h/2)))*C);
   end 
% extend the matrix
   Y = interp1(t_list,Y2(1,:),t_max,'linear','extrap');
   M(jj,:) = Y;
   
end
% calculate the error L2 norm
for kk = 2:length(h_list)
   error_list2(1,kk-1) = norm(M(kk,:)-M(1,:));
end

 % Heun's method 
for jj = 1:length(h_list)
   h = h_list(jj);
   t_list = 0:h:15;
   Y3 = zeros(6,length(t_list));
   Y3(:,1) = [1 0.5 0.25 0 0 -1]';
    
   for ii = 2:length(t_list)
       t1 = t_list(ii);
       t2 = t_list(ii-1);
       Y3(:,ii) = inv(eye(6)-0.75*h*inv(A)*(-B))*(Y3(:,ii-1)+ ...,
           0.25*h*inv(A)*(-B)*Y3(:,ii-1)+0.25*h*inv(A)*(3*sin(2*pi*t1)+sin(2*pi*t2))*C);
   end  
% extend the matrix
  
   Y = interp1(t_list,Y3(1,:),t_max,'linear','extrap');
   M(jj,:) = Y;
    
end
% calculate the error L2 norm
for kk = 2:length(h_list)
   error_list3(1,kk-1) = norm(M(kk,:)-M(1,:));
end
  
loglog(h_list,error_list1,'g','linewidth',2);
hold on
loglog(h_list,error_list2,'b','linewidth',2);
hold on
loglog(h_list,error_list3,'r','linewidth',2);
grid on
hold on

title('Three Second Order Implicit Method');
xlabel('Time step');
ylabel('L_2 norm of error');
legend('trapezoidal rule','midpoint method','Heuns method');