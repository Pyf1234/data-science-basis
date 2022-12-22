%�˴���eΪ���Φ�,���޸Ĵ˴���ֵ����һ���޸ĺ��������ֵ��Ŀǰ����֪�����������ͬ���ĸ�
e=20;
t=0:0.001:30;
n=length(t);
%���庯��������������
fa=@(x,y) -e.*(y.^2-1).*x-y;
fb=@(x) x;
%ŷ�������
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
%������������
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
%���ú���
%ʹ��matlab���õ�ode45ָ��ֱ�����
[t1,y] = ode45(@vdp1,[0:1e-3:30],[0.1; -1]);
%�˴��Ļ�ͼ����д�˺ü�������Ҫ���ڦ�ȡĳЩֵ��ʱ�����ַ����޷�����һ��ͼ����
%��������ʦ���Դ��룬�������а�����һ����ע��ȥ�������ǰ�����һ����ע��
%plot(t,y_e);
%plot(t,y_r,t1,y(:,1))
%plot(t,y_r);
plot(t1,y(:,1));
%plot(t,y_e,t,y_r,t1,y(:,1))
%�˴���y��1����Ϊy��y��2����Ϊy��һ�׵�
function dydt = vdp1(t1,y)
dydt = [y(2); 
		20*(1-y(1)^2)*y(2)-y(1)];%���Ų���ʡ�����Ǳ�����y��2x1������
end

