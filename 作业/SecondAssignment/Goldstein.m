function [alpha, xk, f, k] = Goldstein(fun, grid, x0, dk)
	%
	% Function [alpha, xk, fx, k] = Goldstein(fun, grid, x0, dk)
	% 求出函数fun在x0处以dk为下降方向时的步长alpha，同时返回相对应的下
	% 一个下降点xk以及xk处的函数值fx，k为迭代次数
	% -----------------------------------------------------------
	% 输入: 
	% 	fun 	函数名称(字符变量）
	%	grid 	梯度函数名称(字符变量)
	%	x0		迭代点(列向量)
	%	dk		函数在迭代点处的下降方向(列向量)
	%
	% 输出:
	%	alpha	函数在x0处以dk为下降方向时的下降步长
	%	xk		函数在x0处以dk为下降方向，以alpha为步长
	%			求得的下降点
	%	f	    函数在下降点xk处的函数值
	%	k		求步长算法迭代次数
	% -----------------------------------------------------------
	%
	c = 0.3; 	% 泰勒展开式补足系数，0 < c < 1/2
	alpha = 1; 	% 初始步长为 1
	k = 0; 		% 统计迭代次数
    a = 0; b = inf; % 二分法确定 alpha 值
	gk = feval(grid, x0);	% x0处的梯度值
	fk = feval(fun, x0 + alpha * dk); 	% 函数在下一个迭代点处的目标函数值
	l1 = feval(fun, x0) + c * alpha * gk' * dk; 	% Armjio准则
    l2 = feval(fun, x0) + (1 - c) * alpha * gk' * dk; 	% Armjio准则的补全
	while true
	    if fk > l1
            k = k + 1;
            b = alpha;
            alpha = (a + b) / 2;
            fk = feval(fun, x0 + alpha * dk);
            l1 = feval(fun, x0) + c * alpha * gk' * dk;
            l2 = feval(fun, x0) + (1 - c) * alpha * gk' * dk;
            continue;
        end
        if fk < l2
            k = k + 1;
            a = alpha;
            alpha = min([2 * alpha, (a + b) / 2]);
            fk = feval(fun, x0 + alpha * dk);
            l1 = feval(fun, x0) + c * alpha * gk' * dk;
            l2 = feval(fun, x0) + (1 - c) * alpha * gk' * dk;
            continue;
        end
        break;
	end
	xk = x0 + alpha * dk;	% 下降点
	f = feval(fun, xk);	    % 下降点处函数值
end
