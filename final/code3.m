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
%�������x_r=0��������X��ֵ
t1 = 0:0.01:2;
[t,X] = ode45(@(t,X)vdp1(t,X,A,C),t1,X_0);
X=X';
%��������������ֵ�ֽ�ķ��������ɷַ���
%���ȶ�ԭʼ���ݽ��б�׼��
[m,n]=size(X);
mn=mean(X,2);
X=X-repmat(mn,1,n);
%Ȼ�����Э������󣬼���������ֵ����������
cx=(1/(n-1))*X*X';
[V,D]=eig(cx);
lambda=diag(D);
[dummy,m_arrange]=sort(-1*lambda);
lambda=lambda(m_arrange);
%�������ǽ���ȡ���䷽����������������Ӧ�ı任
V=V(:,1:2);
%��ԭʼ�������任,��͵õ���������Ҫ���������ɷ�z1,z2
Z=V'*X;
figure
plot(t, Z(1,:),'b',t,Z(2,:),'r')
xlabel('t'),ylabel('x')
legend('Z_1','Z_2');
%���Ƿ����ڿ�ʼ�ļ����㣬�����˵����������ĵ㣬��˰�����ҵҪ������˵���ѿ�ͷ�ĸ���ȥ��
Z=Z(:,5:length(t));
t=t(5:length(t));
figure
plot(t, Z(1,:),'b',t,Z(2,:),'r')
xlabel('t'),ylabel('z')
legend('Z_1','Z_2');
%������㵼��
dt=0.01;
dz1_dt=zeros(1,length(t));
dz2_dt=zeros(1,length(t));
n=length(t);
for i=2:n-1
    dz1_dt(i)=(Z(1,i+1)-Z(1,i-1))/(2*dt);
    dz2_dt(i)=(Z(2,i+1)-Z(2,i-1))/(2*dt);
end
dz1_dt(1)=(-3*Z(1,1)+4*Z(1,2)-Z(1,3))/(2*dt);
dz2_dt(1)=(-3*Z(2,1)+4*Z(2,2)-Z(2,3))/(2*dt);
dz1_dt(length(t))=(3*Z(1,length(t))-4*Z(1,length(t)-1)+Z(1,length(t)-2))/(2*dt);
dz2_dt(length(t))=(3*Z(2,length(t))-4*Z(2,length(t)-1)+Z(2,length(t)-2))/(2*dt);
figure
plot(t, dz1_dt,'b',t,dz2_dt,'r')
xlabel('t'),ylabel('x')
legend('dZ1/dt','dZ1/dt');
%������Ŀ������ģ�ͣ�����ֱ���󵼣��������parameter�ı��ʽ��Ȼ�����
sum_z1_z_1=dz1_dt*(Z(1,:)');
sum_z2_z_2=dz2_dt*(Z(2,:)');
sum_z1_z_2=dz2_dt*(Z(1,:)');
sum_z1_z2=Z(1,:)*Z(2,:)';
sum_z2_z_1=dz1_dt*(Z(2,:)');
sum_z2_z2=Z(2,:)*Z(2,:)';
sum_z1_z1=Z(1,:)*Z(1,:)';
A=zeros(2,2);
A(1,1)=(sum_z1_z_1*sum_z2_z2-sum_z2_z_1*sum_z1_z2)/(sum_z1_z1*sum_z2_z2-sum_z1_z2^2);
A(1,2)=(sum_z1_z1*sum_z2_z_1-sum_z1_z_1*sum_z1_z2)/(sum_z1_z1*sum_z2_z2-sum_z1_z2^2);
A(2,1)=(sum_z1_z_2*sum_z2_z2-sum_z2_z_2*sum_z1_z2)/(sum_z1_z1*sum_z2_z2-sum_z1_z2^2);
A(2,2)=(sum_z1_z1*sum_z2_z_2-sum_z1_z_2*sum_z1_z2)/(sum_z1_z1*sum_z2_z2-sum_z1_z2^2);
dz_dt_fit=A*Z;
figure
plot(t, dz1_dt,'b',t,dz_dt_fit(1,:),'r')
xlabel('t'),ylabel('x')
legend('dZ1/dt','dZ1/dt_fit');
figure
plot(t, dz2_dt,'b',t,dz_dt_fit(2,:),'r')
xlabel('t'),ylabel('x')
legend('dZ2/dt','dZ2/dt_fit');


function dxdt = vdp1(t,x,A,C)
dxdt = A*x+C;
end