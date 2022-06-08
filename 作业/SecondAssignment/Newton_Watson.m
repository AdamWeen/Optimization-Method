n = 3; %维数
method = "DFP"
x0 = zeros(n,1); %初始点
H0 = eye(n);    %初始矩阵
x = x0;
H = H0;
X = sym(zeros(n,1));
e = 1e-5;  %精度
step = 0;   %迭代步数
%下面计算目标函数f及其梯度和Hessian矩阵
for i = 1:n
    X(i,1) = sym(['x_',num2str(i)]);
end
f = 0; %函数值
for i = 1:29
    r = 0;
    sum = 0;
    for j = 1:n
        r = r + (j-1)*X(j)*(i/29)^(j-2);
        sum = sum + X(j)*(i/29)^(j-1);
    end
    r = r-sum^2-1;
    f = f + r^2;
end
f = f+X(1)^2+(X(2)-X(1)^2-1)^2;
g = sym(zeros(n,1));
for i = 1:n
    g(i,1) = diff(f,X(i));%函数的偏导数值
end
grad0 = zeros(14,1);%用于在n=2时绘制图像
func = zeros(14,3);
while(1)
    step = step+1;
    g_val_old = eval(subs(g,X,x));
    d = -H*g_val_old;
    alpha = Armijo(f, X, x, g_val_old, d, 0.9, 0.001);
    x_new = x+alpha*d;%按照Armijo准则更新
    g_val_new = eval(subs(g,X,x_new));
    grad = sqrt(g_val_new'*g_val_new);
    if grad <= e
        break;
    end
    s = x_new - x;
    y = g_val_new - g_val_old;
    H = update(H, s, y, method);
    x = x_new;
    disp("迭代代数："+step+"，梯度的模长为"+grad+"，迭代点为("+ ...
        toString(x_new)+"),函数值为"+eval(subs(f,X,x)))
disp("("+toString(x_new)+")")
%     func(step,:) = [x_new' eval(subs(f,X,x))];
%     grad0(step,1) = grad;
end
disp("-------------------------------------------------------------------")
disp("方法为:"+method+";迭代次数："+step+",极小点为("+toString(x_new)+ ...
    "),函数值为"+eval(subs(f,X,x)))
% grad0(step,1) = grad;
% func(step,:) = [x_new' eval(subs(f,X,x))];
% plot(func(:,1),func(:,2),'r')
% hold on
% plot(func(end,1),func(end,2),'r*')
% plot(func(1,1),func(1,2),'ro')
% legend("","ending point","starting point")
% hold off
% plot(grad0,"b")

function H1 = update(H, s, y, method)
%SR1,DFP,BFGS方法
switch method
    case "SR1"
        H1 = H + (s-H*y)*((s-H*y)')/((s-H*y)'*y);
    case "DFP"
        H1 = H + s*s'/(s'*y) - H*y*y'*H/(y'*H*y);
    case "BFGS"
        H1 = H + (1+(y'*H*y)/(y'*s))*((s*s')/(y'*s))-((s*y'*H+H*y*s')/(y'*s));
end
end

function alpha = Armijo(f, X, x0, g0, d, sigma, rho)
%Armijo准则
alpha = 1;
f0 = subs(f,X,x0);
while(1)
    x = x0+alpha*d;
    fk = subs(f,X,x);
    if fk <= f0 + rho*g0'*d*alpha
        break
    end
    alpha = alpha*sigma;
    if alpha <= 10e-2
        break
    end
end
end
function s=toString(x)
[r,~] = size(x);
s = "";
for i = 1:r
    s = s + num2str(x(i));
    if i ~= r
        s = s+",";
    end
end
end
