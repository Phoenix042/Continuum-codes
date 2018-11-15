clc;
%prompt='Enter value of X-high: ';
xh=160;%input(prompt);
%prompt='Enter value of Y-high: ';
yh=160;%input(prompt);
X=[0:1:xh];
Y=[0:1:yh];
mu=1e20;%323.1923e9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% mu and lambda has to be in order of e20 inn order to get correct values
%%% of stresses which can have sufficient low error with cutoff radius 1.5A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=1e20;
E=4*mu*(mu+lambda)/(2*mu+lambda);
nu=lambda/(2*mu+lambda);
k=E*(3.35e-10)^3/(12*(1-nu^2));
k1=k*(1+nu)/E;
xc=xh/2;
yc=yh/2;
s_xx=zeros(length(X),length(Y));
s_yy=zeros(length(X),length(Y));
s_xy=zeros(length(X),length(Y));
f=zeros(length(X),length(Y));
for i=1:length(X)
    for j=1:length(Y)
        x=((xh/2)-X(i))*1e-10;
        y=((yh/2)-Y(j))*1e-10;
        r=sqrt(x^2+y^2);
        if r>(2e-10)
            s_xx(i,j)=(1e-9)*2*mu*k1*(y^2-x^2)/r^4;
            s_yy(i,j)=-(1e-9)*2*mu*k1*(y^2-x^2)/r^4;
            s_xy(i,j)=-(1e-9)*2*mu*x*y*((2*k1/r^2)-(1/6))/r^2;
            f(i,j)=sqrt((x^2+y^2)/3)*1e10;
        else
            s_xx(i,j)=NaN;
            s_yy(i,j)=NaN;
            s_xy(i,j)=NaN;
            f(i,j)=NaN;
        end
    end
end
s1=max(max(s_xx));
s2=max(max(s_yy));
s3=max(max(s_xy));
s4=max(max(f));
s5=min(min(s_xx));
s6=min(min(s_yy));
s7=min(min(s_xy));
s8=min(min(f));
disp('Disclination Stresses using displacement field')
fprintf('s_11 max:%d \t min:%d \n',s1,s5)
fprintf('s_22 max:%d \t min:%d \n',s2,s6)
fprintf('s_12 max:%d \t min:%d \n',s3,s7)
fprintf('Z max:%d \t min:%d \n',s4,s8)
surf(X,Y,s_xx)
colorbar;
view(2)
%print('s_11','-dpdf','-r1200')
hold on;
figure(2)
surf(X,Y,s_yy)
colorbar;
view(2)
%print('s_22','-dpdf','-r1200')
hold on;
figure(3)
surf(X,Y,s_xy)
colorbar;
view(2)
%print('s_12','-dpdf','-r1200')
hold on;
figure(4)
surf(X,Y,f)
colorbar;
view(2)
%print('s_21','-dpdf','-r1200')
hold on;