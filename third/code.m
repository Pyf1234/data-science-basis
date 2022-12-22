%此处的e为超参ε,若修改此处的值，请一并修改函数里面的值，目前还不知道如何让两个同步的改
e=20;
t=0:0.001:30;
n=length(t);
%定义函数，方便后面调用
fa=@(x,y) -e.*(y.^2-1).*x-y;
fb=@(x) x;
%欧拉法求解
x_e=zeros(n);
x_e=x_e(:,1);
y_e=zeros(n);
y_e=y_e(:,1);
x_e(1)=-1;
y_e(1)=0.1;
for j=2:n
    x_e(j)=x_e(j-1)+fa(x_e(j-1),y_e(j-1)).*0.5;
    y_e(j)=y_e(j-1)+fb(x_e(j-1)).*0.5;
end
%龙格库塔法求解
x_r=zeros(n);
x_r=x_r(:,1);
y_r=zeros(n);
y_r=y_r(:,1);
x_r(1)=-1;
y_r(1)=0.1;
for j=2:n
    fa1=fa(x_r(j-1),y_r(j-1));
    fb1=fb(x_r(j-1));
    fa2=fa(x_r(j-1)+0.25.*fa1,y_r(j-1)+0.25.*fb1);
    fb2=fb(x_r(j-1)+0.25.*fa1);
    fa3=fa(x_r(j-1)+0.25.*fa2,y_r(j-1)+0.25.*fb2);
    fb3=fb(x_r(j-1)+0.25.*fa2);
    fa4=fa(x_r(j-1)+0.5.*fa3,y_r(j-1)+0.5.*fb3);
    fb4=fb(x_r(j-1)+0.5.*fa3);
    x_r(j)=x_r(j-1)+(fa1+2.*fa2+2.*fa3+fa4).*(0.5/6);
    y_r(j)=y_r(j-1)+(fb1+2.*fb2+2.*fb3+fb4).*(0.5/6);
end
%内置函数
%使用matlab内置的ode45指令直接求解
[t1,y] = ode45(@vdp1,[0:1e-3:30],[0.1; -1]);
%此处的画图函数写了好几个，主要是在ε取某些值的时候，三种方法无法画在一个图里面
%若助教老师调试代码，可以自行把其中一部分注释去掉，或是把其中一部分注释
%plot(t,y_e);
%plot(t,y_r,t1,y(:,1))
%plot(t,y_r);
plot(t1,y(:,1));
%plot(t,y_e,t,y_r,t1,y(:,1))
%此处的y（1）即为y，y（2）即为y的一阶导
function dydt = vdp1(t1,y)
dydt = [y(2); 
		20*(1-y(1)^2)*y(2)-y(1)];%括号不可省，不是变量，y是2x1列向量
end


