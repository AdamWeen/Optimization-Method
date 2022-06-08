function [output] = func_618(f,a0,b0,epsilon)
%用0.618法计算函数极小值的函数
%   提供函数形式、区间左右端点以及误差
lambda = [];mu = [];a = [];b=[];
a=[a,a0];b=[b,b0];k=1;
lambda =  [lambda,a0+0.382*(b0-a0)] ;mu(end+1) = a0+0.618*(b0-a0);
%下一行为画图考虑。Sign用以标记输出哪个误差
kline=zeros(1,100,'int8');error = zeros(1,100,'double');Sign = 0;
while(b(k)-lambda(k)>epsilon&&mu(k)-a(k)>epsilon)
    if(f(lambda(k))>f(mu(k)))
        Sign = 1;
        if(b(k)-lambda(k)<=epsilon)
            output = mu(k);
            break;
        else
        a(k+1)=lambda(k);b(k+1) = b(k);lambda(k+1) = mu(k);
        mu(k+1) = a(k+1)+0.618*(b(k+1)-a(k+1));
        end
    else
        Sign = 2;
        if(mu(k)-a(k)<=epsilon)
            output = lambda(k);
            break;
        else
            a(k+1) = a(k);b(k+1) = mu(k);mu(k+1) = lambda(k);
            lambda(k+1) = a(k+1)+0.382*(b(k+1)-a(k+1));
        end
    end
    kline(k) = k;
    if(Sign==1)
        error(k) = b(k)-lambda(k);
    elseif(Sign==2)
        error(k) = mu(k)-a(k);
    end
    k=k+1;
end
%输出结论
if(b(k)-lambda(k)<=epsilon)
    output = mu(k);
elseif(mu(k)-a(k)<=epsilon)
    output = lambda(k);
end
%加入最后一次迭代的信息
kline(k) = k;
if(Sign==1)
    error(k) = b(k)-lambda(k);
elseif(Sign==2)
    error(k) = mu(k)-a(k);
end
%删除0值
kline(kline==0)=[];
error(error==0)=[];
%绘图
plot(kline,error,'-or');
grid on;
title("Distance between selected points-0.618");
xlabel('Iterations');
ylabel('Distance');

