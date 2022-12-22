R1=20;
R2=10;
R3=25;
R4=10;
R5=30;
R6=40;
V2=0;
V3=200;
A=[R1+R6+R2 -R1 -R2;-R1 R1+R3+R4 -R4;-R2 -R4 R2+R4+R5];
iteration=zeros(2,51);%保存每个V1对应的迭代次数
ave_iteration=zeros(1,2);%保存两种迭代方法的平均迭代次数
result1=zeros(51,3);%保存雅可比迭代的结果
result2=zeros(51,3);%保存高斯赛德尔迭代的结果
times=1;
%下面是雅可比迭代的方法,结果保存在result1中
for V1=0:2:100
   I=[0 0 0];
   j=0;
   while 1
       %这个temp用来存运算之前的值
       temp=I;
       I(1)=(V1-A(1,2)*temp(2)-A(1,3)*temp(3))/A(1,1);
       I(2)=(V2-A(2,1)*temp(1)-A(2,3)*temp(3))/A(2,2);
       I(3)=(V3-A(3,1)*temp(1)-A(3,2)*temp(2))/A(3,3);
       d=sqrt(sum((I-temp).^2));
       if d<1E-6
           break
       end  
       j=j+1;
   end
   iteration(1,times)=j;
   result1(times,:)=temp;
   times=times+1;
end
%下面是高斯赛德尔迭代法，结果保存在result2中
times=1;
for V1=0:2:100
   I=[0 0 0];
   j=0;
   while 1
       %这个temp用来存运算之前的值
       temp=I;
       I(1)=(V1-A(1,2)*temp(2)-A(1,3)*temp(3))/A(1,1);
       I(2)=(V2-A(2,1)*I(1)-A(2,3)*temp(3))/A(2,2);
       I(3)=(V3-A(3,1)*I(1)-A(3,2)*I(2))/A(3,3);
       d=sqrt(sum((I-temp).^2));
       if d<1E-6
           break
       end  
       j=j+1;
   end
   iteration(2,times)=j;
   result2(times,:)=temp;
   times=times+1;
end
%下面计算平均迭代次数
for k=1:51
    ave_iteration(1)=ave_iteration(1)+iteration(1,k);
    ave_iteration(2)=ave_iteration(2)+iteration(2,k);
end
ave_iteration=ave_iteration/51;