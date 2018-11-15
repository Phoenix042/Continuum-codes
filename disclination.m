clear;
clc;
prompt='Enter value of raduis of sheet: ';
r=input(prompt);
k=1;
E=2000;
nu=0.3;
mu=11500000;
count=0;
x=zeros(1+r+((r+1)*(2*3.1/0.1)),1);
y=zeros(1+r+((r+1)*(2*3.1/0.1)),1);
s1=zeros(1+r+((r+1)*(2*3.1/0.1)),1);
s2=zeros(1+r+((r+1)*(2*3.1/0.1)),1);
s3=zeros(1+r+((r+1)*(2*3.1/0.1)),1);
for i=0:0.5:r
    for j=0:0.1:(2*pi)
        count=count+1;
        x(count,1)=i*cos(j);
        y(count,1)=i*sin(j);
        s1(count,1)=-2*mu*k*(1+nu)*cos(2*j)/(i^2*E);
        s2(count,1)=2*mu*k*(1+nu)*cos(2*j)/(i^2*E);
        s3(count,1)=-2*mu*k*(1+nu)*sin(2*j)/(i^2*E);
    end
end
scatter3(x,y,s1,[],s1,'filled');
view(2)
hold on;
colorbar;
figure(2)
scatter3(x,y,s2,[],s2,'filled');
view(2)
hold on;
colorbar;
figure(3)
scatter3(x,y,s3,[],s3,'filled');
view(2)
hold on;
colorbar;