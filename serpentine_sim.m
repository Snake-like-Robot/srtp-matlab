%% author Lzh
close all;
clear;
%% parameters and matrix define
serpentine_para
%% the main loop
delta_t = 0.05;
t_max = 10;
% for c = 0:pi/10:pi
y = zeros(6,1);
% for t = 0:delta_t:t_max
%     t % 打印当前时间，监测程序运行
%     c % 打印当前的c
%     % 更新矩阵 d代表data
%     theta_d = eval(subs(theta));
%     phi_d = eval(subs(phi));
%     phi_dot_d = eval(subs(phi_dot));
%     omegath_d = eval(subs(omegath));
%     J_t_d = eval(subs(J_t));
%     C_t_d = eval(subs(C_t));
%     L_t_d = eval(subs(L_t));
%     B_t_d = D*J_t_d^-1*D';
%     K_t_d = J_t_d^-1*D'*B_t_d^-1;
%     rho = 1/(e'*J_t_d*e);
%     e_rho = rho*e;
%     % 求解R_t S_t Q_t
%     temp_mat = [L_t_d';E']*omegath_d*Df*omegath_d'*[L_t_d E];
%     R_t_d = Dtau+temp_mat(1:n,1:n);
%     S_t_d = temp_mat(1:n,n+1:end);
%     Q_t_d = temp_mat(n+1:end,n+1:end);
%     % 定义微分方程组
%     mI = eye(2)*m;
%     mat1 = e_rho'*R_t_d*e_rho;
%     mat2 = e_rho'*S_t_d;
%     mat3 = S_t_d'*e_rho;
%     mat4 = Q_t_d;
%     mat5 = [e_rho'*R_t_d;S_t_d']*K_t_d*phi_dot_d;
%     [t_period,y_sol] = ode45(@(t,y)odefun(t,y,rho,mI,mat1,mat2,mat3,mat4,mat5),[t t+delta_t],y);
%     plot(y_sol(:,5),y_sol(:,6));
%     hold on;
%     axis equal;
%     drawnow;
%     y = y_sol(end,:)';
% end
% % plot(y_sol(:,5),y_sol(:,6));
% % axis equal;
% % end
% 
% function y_dot = odefun(t,y,rho,mI,mat1,mat2,mat3,mat4,mat5)
%     y_dot(1,:) = -(mat5(1)+mat1*y(1)+mat2*y(2:3))/rho;
%     temp = -mI^-1*(mat5(2:3)+mat3*y(1)+mat4*y(2:3));
%     y_dot(2,:) = temp(1);
%     y_dot(3,:) = temp(2);
%     y_dot(4,:) = y(1);
%     y_dot(5,:) = y(2);
%     y_dot(6,:) = y(3);
% end
%% 绘制增加时间项以后的插值点
% 测试用的数据
x1 = [];
y1 = [];
id = 1;
syms k real;
fx = 1/n*cos(a*cos(k*b/n)+k*c/n);
fy = 1/n*sin(a*cos(k*b/n)+k*c/n);
% for t = 0:delta_t:t_max
%     t % 打印当前时间
% %     Sth_d = eval(subs(Sth));
% %     Cth_d = eval(subs(Cth));
% %     J_t_d = eval(subs(J_t));
% %     C_t_d = eval(subs(C_t));
% %     L_t_d = eval(subs(L_t));
%     x1(id) = symsum(fx,k,1,1+n*omega*t);
%     y1(id) = symsum(fy,k,1,1+n*omega*t);
%     id = id + 1;
% end
% plot(x1,y1,'o');
% hold on;
% x2 = [];
% y2 = [];
% id = 1;
% for t = 0:delta_t:t_max
%     t % 打印当前时间
% %     Sth_d = eval(subs(Sth));
% %     Cth_d = eval(subs(Cth));
% %     J_t_d = eval(subs(J_t));
% %     C_t_d = eval(subs(C_t));
% %     L_t_d = eval(subs(L_t));
%     x2(id) = symsum(fx,k,1,1+n*omega*t+10);
%     y2(id) = symsum(fy,k,1,1+n*omega*t+10);
%     id = id + 1;
% end
% plot(x2,y2,'*');
% hold on;

syms k real;
fx = 1/n*cos(a*cos(k*b/n)+k*c/n);
fy = 1/n*sin(a*cos(k*b/n)+k*c/n);
for i = 1:(n+1)*2
    x1(i) = symsum(fx,k,1,i);
    y1(i) = symsum(fy,k,1,i);
end
plot(x1,y1,x1,y1,'o');
axis equal;
hold on;
L = sqrt((y1(2)-y1(1))^2+(x1(2)-x1(1))^2);
for offset = 0:0

plot(x1(1+offset:n+1+offset),y1(1+offset:n+1+offset),'LineWidth',1.5);
hold on;

xc = 0;yc = 0;theta_eva = 0;
for i = 1+offset:n+offset
    xc = xc + 0.5*(x1(i)+x1(i+1))/n;
    yc = yc + 0.5*(y1(i)+y1(i+1))/n;
%     temp1 = (y1(i+1)-y1(i))/(x1(i+1)-x1(i));
%     temp2 = (y1(2+offset)-y1(1+offset))/(x1(2+offset)-x1(1+offset));
%     theta_eva = theta_eva + 1/n*(atan2(abs(temp1),sign(temp1))-atan2(abs(temp2),sign(temp2)));
    theta_eva = theta_eva + 1/n*(my_atan(x1(i),y1(i),x1(i+1),y1(i+1))-my_atan(x1(1+offset),y1(1+offset),x1(2+offset),y1(2+offset)));
end
theta_eva = theta_eva + my_atan(x1(1+offset),y1(1+offset),x1(2+offset),y1(2+offset));

%点到直线的距离
q1 = [xc yc];
q2 = [xc+cos(theta_eva) yc+sin(theta_eva)];
p = [x1(n+1+offset) y1(n+1+offset)];
d = abs(det([q2-q1;p-q1]))/norm(q2-q1); 

% 判断点在直线上方还是下方
if tan(theta_eva)*(x1(n+1+offset)-xc)+yc <= y1(n+1+offset)%点在直线上方
    
end

% plot(xc,yc,'*');
% hold on;
% plot([xc,xc+0.7*cos(theta_eva)],[yc,yc+0.7*sin(theta_eva)],'LineWidth',1.0);
% hold on;

end