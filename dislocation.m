%clear all;
clc;
b=2.442191e-10;
c=1e-5;
% prompt='Enter value of X-high: ';

xh=800;%input(prompt);
% prompt='Enter value of Y-high: ';
yh=800;%input(prompt);
X=0:10:xh;
Y=0:10:yh;
%E=840.3e9;
nu=0.12;
mu=475.77e9;%323.1923e9;%384.61538e9;
gamma=6*mu;
%lambda=576.92307e9;
l_1=sqrt(c/(2*mu*(1+nu)));
l_2=sqrt(c*(mu+gamma)/(4*mu*gamma));
A=(mu*(1+nu)*b)/(2*pi);
B=(mu*gamma*b)/(pi*(mu+gamma));
s_xx=zeros(length(X),length(Y));
s_yy=zeros(length(X),length(Y));
s_xy=zeros(length(X),length(Y));
s_yx=zeros(length(X),length(Y));
for i=1:length(X)
    for j=1:length(Y)
        x=((xh/2)-X(i))*1e-10;
        y=((yh/2)-Y(j))*1e-10;
        r=sqrt(x^2+y^2);
        p1=y^2+3*(x^2);
        p2=4*(l_1^2)/r^2;
        p3=y^2-3*(x^2);
        p4=2*(y^2)*r/l_1;
        p5=2*(x^2)*r/l_1;
        p6=-p3*4*(l_2^2)/(r^2);
        p7=2*(x^2)*r/l_2;
        p8=4*(l_2^2)/(r^2);
        p9=(x^2-3*(y^2));
        p10=(x/r^4);
        p11=(y/r^4);
        if r>(20e-10)
            s_xx(i,j)=-p11*(A*(p1+p2*p3-p4*bessely(1,(r/l_1))-p3*2*bessely(2,(r/l_1)))-B*(x^2-y^2-p6+p7*bessely(1,(r/l_2))-p3*2*bessely(2,(r/l_2))))/1e9;
            s_yy(i,j)=-p11*(A*(y^2-x^2-p2*p3-p5*bessely(1,(r/l_1))+p3*2*bessely(2,(r/l_1)))+B*(x^2-y^2-p6+p7*bessely(1,(r/l_2))-p3*2*bessely(2,(r/l_2))))/1e9;
            s_xy(i,j)=p10*(A*(x^2-y^2-p2*p9-(2*y^2*r*bessely(1,(r/l_1))/l_1)+(2*p9*bessely(2,(r/l_1))))+B*(x^2+3*(y^2)+p8*p9-(2*x^2*r*bessely(1,(r/l_2))/l_2)-(2*p9*bessely(2,(r/l_2)))))/1e9;
            s_yx(i,j)=p10*(A*(x^2-y^2-p2*p9-(2*y^2*r*bessely(1,(r/l_1))/l_1)+(2*p9*bessely(2,(r/l_1))))-B*(x^2-y^2-p8*p9-(2*y^2*r*bessely(1,(r/l_2))/l_2)+(2*p9*bessely(2,(r/l_2)))))/1e9;;
        else
            s_xx(i,j)=NaN;
            s_yy(i,j)=NaN;
            s_xy(i,j)=NaN;
            s_yx(i,j)=NaN;
        end
    end
end
% surf(X,Y,s_xx)
% colorbar;
% title('\sigma_{11} ','Color','r','FontSize',16)
% xlabel('Position in X direction','FontSize',12,'Color','b')
% ylabel('Position in Y direction','FontSize',12,'Color','b')
% zlabel('\sigma_{11}','FontSize',12,'Color','b')
% view(2)
% print('s_11','-dpdf','-r1200')
% hold on;
% figure(2)
% surf(X,Y,s_yy)
% colorbar;
% title('\sigma_{22} ','Color','r','FontSize',16)
% xlabel('Position in X direction','FontSize',12,'Color','b')
% ylabel('Position in Y direction','FontSize',12,'Color','b')
% zlabel('\sigma_{22}','FontSize',12,'Color','b')
% view(2)
% print('s_22','-dpdf','-r1200')
% hold on;
% figure(3)
% surf(X,Y,s_xy)
% colorbar;
% title('\sigma_{12} ','Color','r','FontSize',16)
% xlabel('Position in X direction','FontSize',12,'Color','b')
% ylabel('Position in Y direction','FontSize',12,'Color','b')
% zlabel('\sigma_{12}','FontSize',12,'Color','b')
% view(2)
% print('s_12','-dpdf','-r1200')
% hold on;
% figure(4)
% surf(X,Y,s_yx)
% colorbar;
% title('\sigma_{21} ','Color','r','FontSize',16)
% xlabel('Position in X direction','FontSize',12,'Color','b')
% ylabel('Position in Y direction','FontSize',12,'Color','b')
% zlabel('\sigma_{21}','FontSize',12,'Color','b')
% view(2)
% print('s_21','-dpdf','-r1200')
% hold on;
s1=max(max(s_xx));
s2=max(max(s_yy));
s3=max(max(s_xy));
s4=max(max(s_yx));
s5=min(min(s_xx));
s6=min(min(s_yy));
s7=min(min(s_xy));
s8=min(min(s_yx));
fprintf('s_11 max:%d \t min:%d \n',s1,s5)
fprintf('s_22 max:%d \t min:%d \n',s2,s6)
fprintf('s_12 max:%d \t min:%d \n',s3,s7)
fprintf('s_21 max:%d \t min:%d \n',s4,s8)