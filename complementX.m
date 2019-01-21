function x=complementX(x,parity)
%作者：赵振华 
%-- Function File: x=complementX(x,parity)
%输出的x将是按照宇称补齐坐标轴左边的部分之后完整的x
%如果parity==1，则x为奇宇称
%如果parity==2，则x为偶宇称
%输出数据
% n=numel(x);
% x1=x(n:-1:1);
% x1=x1(1:n-1);
%修改：
%2012-06-09 自动判断是行向量还是列向量
x1=x(end:-1:1);
x1=x1(1:end-1);
[h,l]=size(x);
if (h>1&&l>1)%判断是否是单纯的向量
    error('complementX:HORL', 'Wrong dimension of X')
end
if parity==1
    %Odd:
    if h>l %列向量
        x = [x1*(-1); x(:,1)];
    else %行向量
        x = [x1*(-1) x(1,:)];
    end
elseif parity==2
    %Even:
    if h>l %列向量
        x = [x1; x(:,1)];
    else %行向量
        x = [x1 x(1,:)];
    end
else
    error('\n The value of EorO is Wrong.\n')
end