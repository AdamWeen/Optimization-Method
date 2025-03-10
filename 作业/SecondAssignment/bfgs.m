function [x,val,k]=bfgs(fun,gfun,x0,varargin)
k=0;
maxk=500;
rho=0.55;
sigma=0.4;
e=1e-5;%精度
n=length(x0);
Hk=eye(n);
while(k<maxk)
    gk=feval(gfun,x0,varargin{:});
    if(norm(gk)<e),break;end
    dk=-Hk\gk;%在后面会更新Hk
    m=0;
    mk=0;
    while(m<20)
        s=feval(fun,x0+rho^m*dk,varargin{:});
        a=feval(fun,x0)+sigma*rho^m*gk'*dk;
        if(s<a)
         mk=m;
         break;
        end
         m=m+1;
        end
        x=x0+dk*rho^mk;
        sk=x-x0;
        yk=feval(gfun,x,varargin{:})-gk;
        if(yk'*sk>0)
          Hk=Hk-(Hk*sk*sk'*Hk)/(sk'*Hk*sk)+(yk*yk')/(yk'*yk);
        end
    x0=x;
    k=k+1;
end
val=feval(fun,x0,varargin{:});
end