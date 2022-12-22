%给参数赋值
m_s=313;
m_tyre=39;
m_hub=7.7;
m_sd=1.2;
m_u=m_tyre+m_hub+m_sd;
k_s=157614;
c_s= 5792;
k_t=2.75e5;
c_t=300;
h_s0=0.7;
h_u0=0.2;
A=[0 1 0 0;-k_s/m_s -c_s/m_s k_s/m_s c_s/m_s;0 0 0 1;k_s/m_u c_s/m_u -(k_s+k_t)/m_u -(c_s+c_t)/m_u];
B=[0 0;0 0;0 0;k_t/m_u c_t/m_u];
C=[0;k_s*(h_s0-h_u0)/m_s;0;(-k_s*(h_s0-h_u0)+k_t*h_u0)/m_u];
%变量初始化
X_0=[0.66;0;0.42;0];
%求得problem1中的x_s的数据
t1 = 0:0.01:2;
[t,X] = ode45(@(t,X)vdp1(t,X,A,C),t1,X_0);
X_s=X(:,1)';
%梯度下降，求最优的c1,c2,c3,这里是用三角函数加一个直流分量拟合，即f1(x)
max_iteration=140000;
learning_rate=0.05;
c=[1;1;1];
for k = 1:max_iteration
    grad=get_grad1(c(1),c(2),c(3),X_s,t1);
    c=c-learning_rate*grad;
    if norm(grad, 2) < 1e-6
        break
    end
end
X_s_=c(1)*sin(c(2)*t1)+c(3);
figure
plot(t, X_s,'b',t,X_s_,'r')
xlabel('t'),ylabel('x')
legend('x_s','x_s(fit)');
%梯度下降，求最优的a1,a2,这里是用直线拟合
max_iteration1=1000;
learning_rate1=0.05;
a=[1;1];
for k1 = 1:max_iteration1
    grad=get_grad2(a(1),a(2),X_s,t1);
    a=a-learning_rate1*grad;
    if norm(grad, 2) < 1e-6
        break
    end
end
X_s_=a(1)*t1+a(2);
figure
plot(t, X_s,'b',t,X_s_,'r')
xlabel('t'),ylabel('x')
legend('x_s','x_s(fit)');
%梯度下降，求最优的b1,b2,b3,这里是用抛物线拟合
max_iteration2=6000;
learning_rate2=0.05;
b=[1;1;1];
for k2 = 1:max_iteration2
    grad=get_grad3(b(1),b(2),b(3),X_s,t1);
    b=b-learning_rate2*grad;
    if norm(grad, 2) < 1e-6
        break
    end
end
X_s_=b(1)*(t1.^2)+b(2)*t1+b(3);
figure
plot(t, X_s,'b',t,X_s_,'r')
xlabel('t'),ylabel('x')
legend('x_s','x_s(fit)');
%内置函数求最优参数c1,c2,c3,c4,c5,c6，这里是用f2(x)拟合
%根据我们的估计,c2=c4=0.07,c1=-3.7,c3=c5=18,c6=0.7
figure
plot(t, X_s)
xlabel('t'),ylabel('x')
legend('x_s');
%上面的原曲线用于确定初始值
c_0=[-3.7,0.07,18,0.07,18,0.7];
c_=fminunc(@(c_) fun(c_,X_s,t1),c_0);
X_s_=exp(c_(1)*t1).*(c_(2)*sin(c_(3)*t1)+c_(4)*cos(c_(5)*t1))+c_(6);
figure
plot(t, X_s,'b',t,X_s_,'r')
xlabel('t'),ylabel('x')
legend('x_s','x_s(fit)');




%这是problem1里的函数，直接搬过来
function dxdt = vdp1(t,x,A,C)
dxdt = A*x+C;
end

%这个函数用来求误差关于三个参数c1,c2,c3的梯度
function grad=get_grad1(c1,c2,c3,x_s,t)
n=length(t);
df_dc1=0;
df_dc2=0;
df_dc3=0;
for i=1:n
    df_dc1=df_dc1+(c1*sin(c2*t(i))+c3-x_s(i))*sin(c2*t(i));
    df_dc2=df_dc2+(c1*sin(c2*t(i))+c3-x_s(i))*c1*c2*cos(c2*t(i))*t(i);
    df_dc3=df_dc3+(c1*sin(c2*t(i))+c3-x_s(i));
end
df_dc1=df_dc1*2/n;
df_dc2=df_dc2*2/n;
df_dc3=df_dc3*2/n;
grad=[df_dc1;df_dc2;df_dc3];
end
%这个函数用来求误差关于a1,a2的梯度
function grad=get_grad2(a1,a2,x_s,t)
n=length(t);
df_da1=0;
df_da2=0;
for i=1:n
    df_da1=df_da1+(a1*t(i)+a2-x_s(i))*t(i);
    df_da2=df_da2+(a1*t(i)+a2-x_s(i));
end
df_da1=df_da1*2/n;
df_da2=df_da2*2/n;
grad=[df_da1;df_da2];
end
%这个函数用来求误差关于三个参数b1,b2,b3的梯度
function grad=get_grad3(b1,b2,b3,x_s,t)
n=length(t);
df_db1=0;
df_db2=0;
df_db3=0;
for i=1:n
    df_db1=df_db1+(b1*(t(i)^2)+b2*t(i)+b3-x_s(i))*(t(i)^2);
    df_db2=df_db2+(b1*(t(i)^2)+b2*t(i)+b3-x_s(i))*t(i);
    df_db3=df_db3+(b1*(t(i)^2)+b2*t(i)+b3-x_s(i));
end
df_db1=df_db1*2/n;
df_db2=df_db2*2/n;
df_db3=df_db3*2/n;
grad=[df_db1;df_db2;df_db3];
end
%这个函数用于计算误差关于c1,c2,c3,c4,c5,c6的函数值
function error=fun(c,x_s,t)
n=length(t);
error=0;
for i=1:n
    error=error+(exp(c(1)*t(i))*(c(2)*sin(c(3)*t(i))+c(4)*cos(c(5)*t(i)))+c(6)-x_s(i))^2;
end
error=error/n;
end