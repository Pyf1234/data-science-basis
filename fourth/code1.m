%ø… ”ªØ
% [x1,x2] = meshgrid(-3:.1:3, -3:.1:3);
% z=(3.65.*(x1.^2)+5.84.*(x2.^2)+2.*0.37.*(x1.*x2))/4+1.2.*x1+4.*x2;
% surf(x1,x2,z);

%exact line search
Q=[3.65 0.37;0.37 5.84];
b=[1.2;4];
x0=[1;1];
x=x0;
f=zeros(1,10);
x_all=zeros(2,10);
max1=0;
for i=1:10
    d_f=0.5*(Q')*x+b;
    if norm(d_f, 2) < 1e-6
        break
    end
    tau=(0.5*(d_f.')*(Q.')*x+(d_f.')*b)/(0.5*(d_f.')*(Q.')*(d_f));
    %tau=0.01;
    x=x-tau.*d_f;
    x_all(:,i)=x;
    f(i)=0.25*(x.')*Q*x+(b.')*x;
    max1=max1+1;
end
i=1:1:10;
i=i(1:max1);
f=f(1:max1);
%plot(i,f);
%backtracking line search
x_=[1;1];
max_num=50;
f_=zeros(1,max_num);
x_all_=zeros(2,max_num);
max2=0;
for i=1:max_num
    d_f=0.5*(Q')*x_+b;
    if norm(d_f, 2) < 1e-6
        break
    end
    step= 1; 
    alpha=0.2; 
    beta =0.5;
    while 1
       x_new= x_- d_f*step;
       f_old=0.25*(x_.')*Q*x_+(b.')*x_;
       f_new=0.25*(x_new.')*Q*x_new+(b.')*x_new;
       if f_new<f_old+ alpha*step*d_f'*(-d_f)
          x_= x_new;
          break
       else
          step=step*beta;
       end
    end
    f_(i)=f_new;
    x_all_(:,i)=x_;
    max2=max2+1;
end
i=1:1:max_num;
i=i(1:max2);
f_=f_(1:max2);
plot(i,f_);
[x_b, f_b] = fminunc(@(x) 0.25 * x' * Q * x + b' * x, x0);
