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
a(1)=(-3*v(1)+4*v(2)-v(3))/(2*dt);
a(n)=(3*v(n)-4*v(n-1)+v(n-2))/(2*dt);
for j=2:n-1
    jerk(j)=(v(j+1)+v(j-1)-2.*v(j))/(dt.^2);
end
jerk(1)=(-v(4)+4.*v(3)-5.*v(2)+2.*v(1))/(dt.^2);
jerk(n)=(-v(n-3)+4.*v(n-2)-5.*v(n-1)+2.*v(n))/(dt.^2);
figure(1),plot(t,v,t,a,t,jerk);