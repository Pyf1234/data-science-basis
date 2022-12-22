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

%1.ode45���΢�ַ�����
%·������Ϊx_r=0
t1 = 0:0.01:2;
[t,X] = ode45(@(t,X)vdp1(t,X,A,C),t1,X_0);
figure
plot(t, X(:,1),'b',t,X(:,3),'r')
xlabel('t'),ylabel('x')
legend('x_s','x_u');
%·������Ϊx_r=0.2sin(0.5t)
t1 = 0:0.01:2;
[t,X] = ode45(@(t,X)vdp2(t,X,A,B,C),t1,X_0);
figure
plot(t, X(:,1),'b',t,X(:,3),'r')
xlabel('t'),ylabel('x')
legend('x_s','x_u');
%·������Ϊ0.01sin(25t)+0.02sin(50t)
t1 = 0:0.01:2;
[t,X] = ode45(@(t,X)vdp3(t,X,A,B,C),t1,X_0);
figure
plot(t, X(:,1),'b',t,X(:,3),'r')
xlabel('t'),ylabel('x')
legend('x_s','x_u');

%2.ǰ��ŷ��������ŷ�����x_r=0ʱ�ķ�����Ľ�
%ǰ��ŷ��
t1 = 0:0.01:2;
X=zeros(4,length(t1));
x=X_0;
i=1;
delta_t=0.01;
for t1 = 0:0.01:2
    dxdt = A*x+C;
    X(:,i)=x;
    x=x+dxdt*delta_t;
    i=i+1;
end
figure
plot(t, X(1,:),'b',t,X(3,:),'r')
xlabel('t'),ylabel('x')
legend('x_s','x_u');
%����ŷ��
t1 = 0:0.01:2;
X=zeros(4,length(t1));
x=X_0;
i=1;
delta_t=0.01;
I=eye(size(A));
a=I-delta_t*A;
b=delta_t*C;
for t1 = 0:0.01:2
    X(:,i)=x;
    x=a\(x+b);
    i=i+1;
end
figure
plot(t, X(1,:),'b',t,X(3,:),'r')
xlabel('t'),ylabel('x')
legend('x_s','x_u');

%3.�����̬��
%�������ſɱȵ���
X(1:2,1)=[X_0(1);X_0(3)];
i=2;
threshold=1e-6;
while true
    X(1,i)=(-C(2)-A(2,3)*X(2,i-1))/A(2,1);
    X(2,i)=(-C(4)-A(4,1)*X(1,i-1))/A(4,3);
    error = norm(X(1:2,i) - X(1:2,i-1), Inf);
    if error<threshold
        break
    end
    i=i+1;
end
t=1:i;
figure
plot(t, X(1,1:i),'b',t,X(2,1:i),'r')
xlabel('iteration'),ylabel('x')
legend('x_s','x_u');
%���ף���˹���¶�����
X(1:2,1)=[X_0(1);X_0(3)];
i=2;
threshold=1e-6;
while true
    X(1,i)=(-C(2)-A(2,3)*X(2,i-1))/A(2,1);
    X(2,i)=(-C(4)-A(4,1)*X(1,i))/A(4,3);
    error = norm(X(1:2,i) - X(1:2,i-1), Inf);
    if error<threshold
        break
    end
    i=i+1;
end
t=1:i;
figure
plot(t, X(1,1:i),'b',t,X(2,1:i),'r')
xlabel('iteration'),ylabel('x')
legend('x_s','x_u');
%����̬�����֤
t1 = 0:0.1:10;
[t,X] = ode45(@(t,X)vdp1(t,X,A,C),t1,X_0);
figure
plot(t, X(:,1),'b',t,X(:,3),'r')
xlabel('t'),ylabel('x')
legend('x_s','x_u');

function dxdt = vdp1(t,x,A,C)
dxdt = A*x+C;
end
function dxdt = vdp2(t,x,A,B,C)
U=[0.2*sin(0.5*t);0.1*cos(0.5*t)];
dxdt = A*x+B*U+C;
end
function dxdt = vdp3(t,x,A,B,C)
U=[0.01*sin(25*t)+0.02*sin(50*t);0.25*cos(25*t)+cos(50*t)];
dxdt = A*x+B*U+C;
end

