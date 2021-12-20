% close all;
clear;
%% robot parameters
robot_para
%% serpentine locomotion
serpentine_fun
%% define matrix
robot_matrix_def
%% the main loop
t = 0;
delta_t = 0.1;
t_max = 10;
y = zeros(6,1);
for t = 0:delta_t:t_max
    t
%% 计算矩阵与数值
% 带花边的J记为J_t，以此类推
    theta_dot = eval(subs(pinv(D)*phi_dot));
    theta = eval(pinv(D)*cellfun(@subs,phi));

    Sth = diag(sin(theta));
    Cth = diag(cos(theta));
    J_t = eval(subs(J+Sth*H*Sth+Cth*H*Cth));
    C_t = eval(subs(Sth*H*Sth-Cth*H*Cth));
    L_t = eval(subs([Sth*N' -Cth*N']'));
    J_t_inv = inv(J_t);
    B_t = D*J_t_inv*D';
    B_t_inv = inv(B_t);
    K_t = J_t_inv*D'*B_t_inv;
    rho = 1/(e'*J_t*e);
    e_rho = rho*e;
    omega_th = [Cth -Sth;Sth Cth];
    
    omega_th_data = eval(subs(omega_th));
    temp_mat1 = [L_t';E']*omega_th_data*Df*omega_th_data'*[L_t E];
    temp_mat2 = [Dtau zeros(n,2);zeros(2,n+2)];
    temp_mat3 = temp_mat1+temp_mat2;
    R_t = temp_mat3(1:n,1:n);
    S_t = temp_mat3(1:n,n+1:end);
    Q_t = temp_mat3(n+1:end,n+1:end);
%% 微分方程组
% 记psi_dot = y1 psi = y2 omegax_dot = y3  omegay_dot = y4 omegax = y5
% omegay = y6
% 记mat1 = e_rho'*R_t*e_rho mat2 = e_rho'*S_t mat3 = S_t'*e_rho
    mat1 = e_rho'*R_t*e_rho;
    mat2 = e_rho'*S_t;
    mat3 = S_t'*e_rho;
    mat4 = [e_rho'*R_t;S_t']*K_t*eval(subs(phi_dot));
    % 大写的I是什么东西 (猜测是2*2的单位阵)
    I = eye(2);
    mI = m*I;
    [t,y] = ode45(@(t,y)odefun(t,y,mat1,mat2,mat3,mat4,rho,mI,Q_t),[t t+delta_t],y);
    plot(y(:,5),y(:,6));
    hold on;
    drawnow
%     plot(t,sqrt(y(:,3).^2+y(:,4).^2));
%     hold on;
%     plot(cos(theta(1))*li+y(:,5),sin(theta(1))*li+y(:,6));
%     hold on;
%     drawnow
%     plot(t,eval(subs(phi{2})));
%     hold on;
%     drawnow
%     plot(t,eval(subs(phi{3})));
%     hold on;
%     drawnow
%     plot(t,eval(subs(phi{4})));
%     hold on;
%     drawnow
%     plot(t,eval(subs(phi{5})));
%     hold on;
%     drawnow
    y = y(end,:)';
end
function dydt = odefun(t,y,m1,m2,m3,m4,rho,mI,Q_t)
%ODEFUN 代换之后的微分方程组
%   暂无
dydt = zeros(4,1);
dydt(2) = y(1);
dydt(5) = y(3);
dydt(6) = y(4);
dydt(1) = -(m4(1)+m2*y(3:4)+m1*y(1))/rho;
temp= -mI\(m4(2:end)+Q_t*y(3:4)+m3*y(1));
dydt(3) = temp(1);
dydt(4) = temp(2);
end