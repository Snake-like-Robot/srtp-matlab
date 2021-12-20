close all;
clear;

plot([0 50],[0 0],'LineStyle',"--",'LineWidth',1);
axis equal;
hold on;
plot([0 0],[-30 30],'LineStyle',"--",'LineWidth',1);
axis equal;
hold on;

n = 7;
L = 7;
list = zeros(4,n+1);%x,y,angle(absolute),angle(relative)
h_list = zeros(1,n);
for i = 1:n
    h_list(i) = plot([list(1,i) list(1,i+1)],[list(2,i) list(2,i+1)]);
    hold on;
end
h_dir = plot([0 0],[0 0],'l)
t_ms = 0;
alpha = 30;
beta = 60;
gamma = 0;
omega = 5;
while true
    t_ms = t_ms+1;
    list(1,1) = 0;
    list(2,1) = 0;
    theta_ave = 0;
    for i = 1:n
        list(3,i) = 0;
        list(4,i) = alpha*sin(omega*t_ms+(i-1)*beta/180*pi)+gamma;
        for j = 1:i
            list(3,i) = list(3,i)+list(4,j);
        end
        list(1,i+1) = list(1,i)+L*cos(list(3,i)/180*pi);
        list(2,i+1) = list(2,i)+L*sin(list(3,i)/180*pi);
        set(h_list(i),'XData',[list(1,i) list(1,i+1)],'YData',[list(2,i) list(2,i+1)]);
    end
    for i = 1:n
        theta_ave = theta_ave+list(4,i);
    end
    theta_ave = theta_ave/n;
    drawnow;
    pause(0.05)
end