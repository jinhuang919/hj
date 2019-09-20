N=input('N=');
%对区域进行网格划分
h1=pi/N;
h2=1/N;
x=[0:h1:pi];
y=[0:h2:1];
[X,Y]=meshgrid(x,y);
%精确解
%v=(9+pi^2)\cos(3*X).*sin(pi*Y)
v=zeros(N+1);
for i=1:N
    for j=1:N
        v(i,j)=(9+pi^2)\cos(3*x(i))*sin(pi*y(j));
    end
end
%v=v';
%构造需要的矩阵及常向量
X=X(2:N,2:N);Y=Y(2:N,2:N);
A=(-1/(h2)^2)*eye(N-1);
C=(-1/h1^2)*eye(N-2);
B=2*(1/(h1^2)+1/(h2^2))*eye(N-1)+diag(diag(C),1)+diag(diag(C),-1);
B(1,1)=1/(h1^2)+2/h2^2;
B(N-1,N-1)=B(1,1);
E=diag(diag(eye(N-2)),1);
F=diag(diag(eye(N-2)),-1);
f=cos(3*X).*sin(pi*Y);
%构造系数矩阵
D = kron(eye(N-1),B)+kron(E,A)+kron(F,A);
%求解方程组
u = zeros(N+1,N+1);f=f';
u(2:N,2:N)=reshape(D\f(:),N-1,N-1);%求解内点处函数值
u(:,1)=u(:,2);
u(:,N+1)=u(:,N);
%u=u';
%画出差分解及误差波动图像
error=u-v;
error=error(2:N,2:N);
%figure(1),mesh(x,y,error);
%xlabel x,ylabel y,zlabel error;
%figure(2),mesh(x,y,u)
%xlabel x,ylabel y,zlabel u
%figure(3),mesh(x,y,v)
%xlabel x,ylabel y,zlabel v
Error=norm(error(:),inf)