R1=20;
R2=10;
R3=25;
R4=10;
R5=30;
R6=40;
V2=0;
V3=200;
A=[R1+R6+R2 -R1 -R2;-R1 R1+R3+R4 -R4;-R2 -R4 R2+R4+R5];
I=zeros(51,3);
i=1;
for V1=0:2:100
    V=[V1 V2 V3].';
    solution=A\V;
    I(i,:)=solution.';
    i=i+1;
end
