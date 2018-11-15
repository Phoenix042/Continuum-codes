%clear all;
clc;
%prompt='Enter the magnitude of burgersvector: ';
b=2.442191;%input(prompt);
b=b*1e-10;
%prompt='Enter value of dislocation modulus: ';
c=1e-7;%input(prompt);
%prompt='Enter value of X-high: ';
xh=80;%input(prompt);
%prompt='Enter value of Y-high: ';
yh=80;%input(prompt);
X=[0:1:xh];
Y=[0:1:yh];
E=1e12;
nu=0.3;
mu=384.61538e9;
gamma=mu/2;
lambda=576.92307e9;
l_1=sqrt(c/(2*mu*1.3));
l_2=sqrt(c*1.5/(2*mu));
A=(mu*1.3*b)/(2*pi);
B=(mu*gamma*b)/(pi*mu*1.5);
s_xx=zeros(length(X),length(Y));
s_yy=zeros(length(X),length(Y));
s_xy=zeros(length(X),length(Y));
s_yx=zeros(length(X),length(Y));
for i=1:length(X)
    for j=1:length(Y)
        x=((xh/2)-X(i))*1e-10;
        y=((yh/2)-Y(j))*1e-10;
        r=sqrt((X(i)-(xh/2))^2+(Y(j)-(yh/2))^2)*1e-10;
        p1=y^2+3*x^2;
        p2=4*l_1^2/r^2;
        p3=y^2-3*x^2;
        p4=2*y^2*r/l_1;
        p5=2*x^2*r/l_1;
        p6=-p3*4*l_2^2/r^2;
        p7=2*x^2*r/l_2;
        p8=4*l_2^2/r^2;
        s_xx(i,j)=-(y/r^4)*(A*(p1+p2*p3-p4*bessely(1,(r/l_1))-p3*2*bessely(2,(r/l_1)))-B*(x^2-y^2-p6+p7*bessely(1,(r/l_2))-p3*2*bessely(2,(r/l_2))))/1e9;
        s_yy(i,j)=-(y/r^4)*(A*(y^2-x^2-p2*p3-p5*bessely(1,(r/l_1))+p3*2*bessely(2,(r/l_1)))+B*(x^2-y^2-p6+p7*bessely(1,(r/l_2))-p3*2*bessely(2,(r/l_2))))/1e9;
        s_xy(i,j)=(x/r^4)*(A*(x^2-y^2-p2*(x^2-3*y^2)-(2*y^2*r*bessely(1,(r/l_1))/l_1)+(2*(x^2-3*y^2)*bessely(2,(r/l_1))))+B*(x^2+3*y^2+p8*(x^2-3*y^2)-(2*x^2*r*bessely(1,(r/l_2))/l_2)-(2*(x^2-3*y^2)*bessely(2,(r/l_1)))))/1e9;
        s_yx(i,j)=(x/r^4)*(A*(x^2-y^2-p2*(x^2-3*y^2)-(2*y^2*r*bessely(1,(r/l_1))/l_1)+(2*(x^2-3*y^2)*bessely(2,(r/l_1))))-B*(x^2-y^2-p8*(x^2-3*y^2)-(2*y^2*r*bessely(1,(r/l_2))/l_2)+(2*(x^2-3*y^2)*bessely(2,(r/l_1)))))/1e9;;
    end
end
surf(X,Y,s_xx)
colorbar;
%view(2)
%print('s_11','-dpdf','-r1200')
hold on;
figure(2)
surf(X,Y,s_yy)
colorbar;
%view(2)
%print('s_22','-dpdf','-r1200')
hold on;
figure(3)
surf(X,Y,s_xy)
colorbar;
%view(2)
%print('s_12','-dpdf','-r1200')
hold on;
figure(4)
surf(X,Y,s_yx)
colorbar;
%view(2)
%print('s_21','-dpdf','-r1200')
hold on;