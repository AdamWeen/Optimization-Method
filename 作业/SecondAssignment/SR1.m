function [x,val,k]=SR1(fun,gfun,x0)
%功能：用秩1拟牛顿法求无约束问题 mini f(x)
%输入：fun,gfun分别是目标函数和梯度，x0是初始点
%输出：x,val分别是近似最优点和最优值，k表示迭代次数
k=0;
maxk=500;
rho=0.55;
sigma=0.4;
e=1e-5;%精度
n=length(x0);
Hk=eye(n);
while(k<maxk)
    gk=feval(gfun,x0);
    dk=-Hk*gk;
    if(norm(gk)<=e),break;end
m=0;
mk=0;
while(m<20)
        if(feval(fun,x0+rho^m*dk)<feval(fun,x0)+sigma*rho^m*gk'*dk)
           mk=m;
           break;
        end
       m=m+1;
       
end
x=x0+dk*rho^mk;
sk=x-x0;
yk=feval(gfun,x)-gk;
Hk=Hk+(sk-Hk*yk)*(sk-Hk*yk)'/((sk-Hk*yk)'*yk);
k=k+1;
x0=x;
end
x0=x;
val=feval(fun,x0);
end