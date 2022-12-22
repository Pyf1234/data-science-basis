close all; 
clear all;
%确定区间
dt=0.1;
t=0:dt:4;
n=length(t);
%这里计算速度在每个小区间上的值
v=exp(-0.1.*t).*cos(5.*t) + t.^2 - 0.1.*(t.^4);
for j=2:n-1
    a(j)=(v(j+1)-v(j-1))/(2*dt);
end
%这里计算两个端点值
a(1)=(-3*v(1)+4*v(2)-v(3))/(2*dt);
a(n)=(3*v(n)-4*v(n-1)+v(n-2))/(2*dt);
figure(1),plot(t,v,t,a);

