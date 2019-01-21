function x=complementX(x,parity)
%���ߣ����� 
%-- Function File: x=complementX(x,parity)
%�����x���ǰ�����Ʋ�����������ߵĲ���֮��������x
%���parity==1����xΪ�����
%���parity==2����xΪż���
%�������
% n=numel(x);
% x1=x(n:-1:1);
% x1=x1(1:n-1);
%�޸ģ�
%2012-06-09 �Զ��ж�������������������
x1=x(end:-1:1);
x1=x1(1:end-1);
[h,l]=size(x);
if (h>1&&l>1)%�ж��Ƿ��ǵ���������
    error('complementX:HORL', 'Wrong dimension of X')
end
if parity==1
    %Odd:
    if h>l %������
        x = [x1*(-1); x(:,1)];
    else %������
        x = [x1*(-1) x(1,:)];
    end
elseif parity==2
    %Even:
    if h>l %������
        x = [x1; x(:,1)];
    else %������
        x = [x1 x(1,:)];
    end
else
    error('\n The value of EorO is Wrong.\n')
end