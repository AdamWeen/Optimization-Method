function [output,iteration,outPut] = SR1(x1,epsilon)
%SR1方法(对称秩1修正公式)
%   输入初始迭代点,停机条件,输出搜索结果,迭代步数,每步结果
    m=1000;n=size(x1);B(:,:,1)=eye(n);k=1;outPut=zeors(m,1);
    %索引从1开始,故初值为x1;m是预设迭代次数
    outPut(1)=Watson(x1);x=zeros(n,m);x(1)=x1;
    d=zeros(n,m);g=zeros(n,m);s=zeors(n,m);y=zeros(n,m);
    
    g(:,k)=gradient(Watson,x1);
end