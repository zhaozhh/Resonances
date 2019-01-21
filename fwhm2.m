function [width, flag]= fwhm2(x,y)
%2013.07.07添加
%如果flag==1半高宽存在，如果flag==-1则不存在。
%width =sqrt(right) - sqrt(left);   % 峰的半高宽 %和fwhm的定义不同
%
%
%
%x:坐标
%y: the line of peak
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)
%Modified by Zhen-Hua Zhao
%2010/12/05

[peakValue2,peakPostion]=max(y);
%peakPostion=find(diff(sign(diff(y)))<0)+1;%找出峰的最大值所在的位置。
%if numel(peakPostion)>1
%    error('myApp:fwhm','峰的个数大于一个')
%end
y = y / y(peakPostion);%这样峰的最大值就是1了。
N = length(y);%总点数
lev50 = 0.5;
centerindex=peakPostion;%以第一个峰的最大值的位置为峰的位置
% fprintf('centerindex= %f \n',centerindex)
% fprintf('numel(x)= %f \n',numel(x))


%————————————————————————————————————————
%首先求峰的左边半高宽处的位置
i = centerindex;
if i == 1 %i的最小值为2;
     fprintf('i= %d \n',i)
     error('myApp:fwhm','这不是一个峰！')
end
%fprintf('i= %d \n',i)
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) && (i >2))
    i = i-1;
end    %first crossing is between v(i-1) & v(i)

if i ~= 2 %i的最小值为2，如果不是2说明峰的左边半高宽处的点是存在的
    interp = (lev50-y(i-1)) / (y(i)-y(i-1)); %斜率
    left = x(i-1) + interp*(x(i)-x(i-1));%峰的左边半高宽处的位置
else
    disp('fwhm2: 峰的左边半高宽处的点不存在')
    if nargout==1
        width = Inf;
    else
        width = Inf;
        flag=-1;
    end
    return
end
%————————————————————————————————————————


%————————————————————————————————————————
%进一步求峰的右边半高宽处的位置
i = centerindex+1;    
%fprintf('i=%i',i)%start search for next crossing at center
while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) && (i <= N-1))
    i = i+1;
end
if i ~= N  
   % disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y(i-1)) / (y(i)-y(i-1)); %斜率
    right = x(i-1) + interp*(x(i)-x(i-1));%峰的右边半高宽处的位置
    if nargout==1
        width =sqrt(right) - sqrt(left);   % 峰的半高宽 
    else
        width =sqrt(right) - sqrt(left);
        flag=1;
    end
	if width<0
		error('fwhm2:   width<0  ');
	end
    
else    
    disp('fwhm2: 峰的右边半高宽处的点不存在')
    if nargout==1
        width = Inf;
    else
        width = Inf;
        flag=-1;
    end
   
    return
end

end