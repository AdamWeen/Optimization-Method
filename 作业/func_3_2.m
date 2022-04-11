function [output] = func_3_2(f,a0,b0,epsilon)
%使用三点二次插值函数求函数极小值点
%   输入函数形式、左右端点以及要求误差
%(默认函数在此区间内任意点函数值小于在端点的函数值）
    x1=a0;x2=(a0+b0)/2;x3=b0;f1=f(x1);f2=f(x2);f3=f(x3);
    x=0.5*((x2^2-x3^2)*f1+(x3^2-x1^2)*f2+(x1^2-x2^2)*f3)/((x2-x3)*f1+ ...
    (x3-x1)*f2+(x1-x2)*f3);
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
    if(abs(f(x2)-f(x))<epsilon*f(x2))
        break;
    end
end
if(f(x)<f(2))
    output = x;
else
    output = x2;
end