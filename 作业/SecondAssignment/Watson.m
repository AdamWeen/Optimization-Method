function [outputArg1] = Watson(x,i)
%Watson函数
%   输入列向量x，下标i，输出函数值
n=size(x,1);sum1=0;sum2=0;t=zeros(29,1);
if i<=29&&i>=1
    for k=1:29
        t(i)=k/29;
    end
    for j=2:n
        sum1=sum1+(j-1)*x(j)*t(i)^(j-2);
    end
    for j=1:n
        sum2=sum2+x(j)*t(i)^(j-1);
    end
    outputArg1=sum1-sum2^2-1;
elseif i==30
    outputArg1=x(1);
elseif i==31
    outputArg1=x(2)-x(1)^2-1;
end
end