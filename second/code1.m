close all; 
clear all;
%使用区间左端点函数值
sum=0;
for t=0:0.01:(4-0.01)
    sum=sum+(exp(-0.1*t)*cos(5*t) + (t.^2 - 0.1*(t.^4)))*0.01;
end
%使用区间右端点函数值
sum1=0;
for t=0:0.01:(4-0.01)
    t_=t+0.01;
    sum1=sum1+(exp(-0.1*t_)*cos(5*t_) + (t_.^2 - 0.1*(t_.^4)))*0.01;
end
%使用梯形法则
sum2=0;
for t=0:0.01:(4-0.01)
    t_=t+0.01;
    sum2=sum2+0.5*(exp(-0.1*t_)*cos(5*t_) + (t_.^2 - 0.1*(t_.^4)))*0.01;
    sum2=sum2+0.5*(exp(-0.1*t)*cos(5*t) + (t.^2 - 0.1*(t.^4)))*0.01;
end
%使用内置函数
f=@(x) exp(-0.1.*x).*cos(5.*x)+(x.^2-0.1.*(x.^4));
[sum3,times]=quad(f,0,4);