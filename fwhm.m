function [width, flag]= fwhm(x,y)
% ���flag==1��߿���ڣ����flag==-1�򲻴��ڡ�
% width = x_right - x_left;   % ��İ�߿� 
% x: ����
% y: �����ߵ�ֵ
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)
%Modified by Zhen-Hua Zhao
%2010/12/05

[peakValue2,peakPostion]=max(y);
%peakPostion=find(diff(sign(diff(y)))<0)+1;%�ҳ�������ֵ���ڵ�λ�á�
%if numel(peakPostion)>1
%    error('myApp:fwhm','��ĸ�������һ��')
%end
y = y / y(peakPostion);%����������ֵ����1�ˡ�
N = length(y);%�ܵ���
lev50 = 0.5;
centerindex=peakPostion;%�Ե�һ��������ֵ��λ��Ϊ���λ��
% fprintf('centerindex= %f \n',centerindex)
% fprintf('numel(x)= %f \n',numel(x))


%��������������������������������������������������������������������������������
%����������߰�߿���λ��
i = centerindex;
if i == 1 %i����СֵΪ2;
     fprintf('i= %d \n',i)
     error('myApp:fwhm','�ⲻ��һ���壡')
end
%fprintf('i= %d \n',i)
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) && (i >2))
    i = i-1;
end    %first crossing is between v(i-1) & v(i)

if i ~= 2 %i����СֵΪ2���������2˵�������߰�߿��ĵ��Ǵ��ڵ�
    interp = (lev50-y(i-1)) / (y(i)-y(i-1)); %б��
    left = x(i-1) + interp*(x(i)-x(i-1));%�����߰�߿���λ��
else
    disp('fwhm: �����߰�߿��ĵ㲻����')
    if nargout==1
        width = Inf;
    else
        width = Inf;
        flag=-1;
    end
    return
end
%��������������������������������������������������������������������������������


%��������������������������������������������������������������������������������
%��һ�������ұ߰�߿���λ��
i = centerindex+1;    
%fprintf('i=%i',i)%start search for next crossing at center
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) && (i <= N-1))
    i = i+1;
end
if i ~= N  
   % disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y(i-1)) / (y(i)-y(i-1)); %б��
    right = x(i-1) + interp*(x(i)-x(i-1));%����ұ߰�߿���λ��
    if nargout==1
        width =right - left;   % ��İ�߿� 
    else
        width =right - left;
        flag=1;
    end
	if width<0
		error('fwhm:   width<0  ');
	end
    
else    
    disp('fwhm: ����ұ߰�߿��ĵ㲻����')
    if nargout==1
        width = Inf;
    else
        width = Inf;
        flag=-1;
    end
   
    return
end

end