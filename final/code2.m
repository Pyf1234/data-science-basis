%��������ֵ
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
%������ʼ��
X_0=[0.66;0;0.42;0];
%���problem1�е�x_s������
t1 = 0:0.01:2;
[t,X] = ode45(@(t,X)vdp1(t,X,A,C),t1,X_0);
X_s=X(:,1)';
%�ݶ��½��������ŵ�c1,c2,c3,�����������Ǻ�����һ��ֱ��������ϣ���f1(x)
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
%�ݶ��½��������ŵ�a1,a2,��������ֱ�����
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
%�ݶ��½��������ŵ�b1,b2,b3,�����������������
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
%���ú��������Ų���c1,c2,c3,c4,c5,c6����������f2(x)���
%�������ǵĹ���,c2=c4=0.07,c1=-3.7,c3=c5=18,c6=0.7
figure
plot(t, X_s)
xlabel('t'),ylabel('x')
legend('x_s');
%�����ԭ��������ȷ����ʼֵ
c_0=[-3.7,0.07,18,0.07,18,0.7];
c_=fminunc(@(c_) fun(c_,X_s,t1),c_0);
X_s_=exp(c_(1)*t1).*(c_(2)*sin(c_(3)*t1)+c_(4)*cos(c_(5)*t1))+c_(6);
figure
plot(t, X_s,'b',t,X_s_,'r')
xlabel('t'),ylabel('x')
legend('x_s','x_s(fit)');




%����problem1��ĺ�����ֱ�Ӱ����
function dxdt = vdp1(t,x,A,C)
dxdt = A*x+C;
end

%���������������������������c1,c2,c3���ݶ�
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
%�������������������a1,a2���ݶ�
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
%���������������������������b1,b2,b3���ݶ�
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
%����������ڼ���������c1,c2,c3,c4,c5,c6�ĺ���ֵ
function error=fun(c,x_s,t)
n=length(t);
error=0;
for i=1:n
    error=error+(exp(c(1)*t(i))*(c(2)*sin(c(3)*t(i))+c(4)*cos(c(5)*t(i)))+c(6)-x_s(i))^2;
end
error=error/n;
end