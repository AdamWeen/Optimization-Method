function [output] = func_3_2(f,a0,b0,epsilon)
%使用三点二次插值函数求函数极小值点
%   输入函数形式、左右端点以及要求误差
%(默认函数在此区间内任意点函数值小于在端点的函数值）
    x1=a0;x2=(a0+b0)/2;x3=b0;f1=f(x1);f2=f(x2);f3=f(x3);
    x=0.5*((x2^2-x3^2)*f1+(x3^2-x1^2)*f2+(x1^2-x2^2)*f3)/((x2-x3)*f1+ ...
    (x3-x1)*f2+(x1-x2)*f3);
    %以下为画图考虑
    k=1;kline=zeros(1,100,'int8');error = zeros(1,100,'double');
while(abs(f(x2)-f(x))>=epsilon*f(x2))
    fx=f(x);f1=f(x1);f2=f(x2);f3=f(x3);
    if(x>x2)
        if(fx<=f2)
            x1=x2;x2=x;
        else
            x3=x;
        end
    else
        if(fx<=f2)
            x3=x2;x2=x;
        else
            x1=x;
        end
    end
    x=0.5*((x2^2-x3^2)*f1+(x3^2-x1^2)*f2+(x1^2-x2^2)*f3)/((x2-x3)*f1+ ...
    (x3-x1)*f2+(x1-x2)*f3);
    kline(k)=k;error(k) = abs(f(x2)-f(x));
    if(abs(f(x2)-f(x))<epsilon*f(x2))
        break;
    end
    k=k+1;
end
if(f(x)<f(2))
    output = x;
else
    output = x2;
end
kline(k)=k;error(k) = abs(f(x2)-f(x));
%删除0值
kline(kline==0)=[];
error(error==0)=[];
plot(kline,error,'-or');
grid on;
title("Distance between selected points-三点二次");
xlabel('Iterations');
ylabel('Distance');