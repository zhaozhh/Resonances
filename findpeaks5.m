function [m,m2,ratio,m2mDeltamTau]=findpeaks5(EorO,U,fraction)
%求所有的峰值及最终的谱函数
%2014.01.19,优化了一下，提供给fwhm函数的是一个完整的峰。
%2014.01.16,搜索区间的起始值改为势能的边界值。
%2013.08.13,修改，改正了没找到峰会报错停止运行的问题。
%2013.03.19,改进。利用每个真实的峰都有半高宽的特性。
%2010年12月为函数ferimionSpectrum,Zhen-Hua Zhao, zhaozhh02@gmail.com。
%EorO:    确定方程解得奇偶性，1：奇函数，2：偶函数。
%chiral： 确定费米子的手征性，1，左手费米子，2右手费米子。
%fraction:最小区间和最大区间的比值。
%mm2Deltam:包含对应峰的位置的m m^2和峰的半高宽deltam的矩阵
%ratio:小区间的波函数平方的积分除以大区间的波函数平方的积分
%m：   质量本征值
%m2：  质量的平方（也就是薛定谔方程的本征值）
%npeaks=5;%峰的个数的估计值。
%******************************************************
%在全m^2的区间内找出峰的个数
%******************************************************
global addition


if EorO==1   
    fprintf('\n 奇函数部分\n ');
elseif EorO==2
    fprintf('\n 偶函数部分\n ');
else
    error(' Wrong with EorO\n ')    
end
fprintf('------------------\n');
fprintf('The maximum of the potential: %f \n', max(U(:,2)));
addition=max(U(:,2))*0.5;%增量，m^2的最大值要比是函数的最高值大addition。
fprintf('\n m^2的最大值要比势能的最大值大%.3f   \n ',addition);
fprintf('m的最大值不超过   ：%f \n', sqrt(max(U(:,2))+addition));
fprintf(' \n');
num=input('输入初始扫描的间隔数(默认为500):   ');
if isempty(num)
    num =500;
end

fprintf('初始扫描的间隔数为 :%i \n', num)
limitm2=max(U(:,2))+addition;
first_m2Space =linspace(U(end,2),limitm2,num);
firstRatio=spectrum(EorO,U,first_m2Space,num,fraction);
pointPeakes=find(diff(sign(diff(firstRatio)))<0)+1;%find out all the maximum.
fprintf('可能的峰的个数为  ：%i\n', length(pointPeakes));
disp('大约位置在：');
fprintf('m^2->：%f\n', first_m2Space(pointPeakes));
disp('------------------');

pointsmin=find(diff(sign(diff(firstRatio)))>0)+1;%find out all the minimum.

%fprintf('pointsmin    ：%f \n', numel(pointsmin));

noumberOfPeakes=numel(pointPeakes);
noumberOfValleies=numel(pointsmin);


if noumberOfPeakes>=1 %Is there a maximum?
    if noumberOfPeakes==1
        if noumberOfValleies~=0
            if pointPeakes(1)<pointsmin(1)
                m2begin=U(end,2);
                m2end  =first_m2Space(pointsmin(1));
                
            else
                m2begin=first_m2Space(pointsmin(1));
                m2end  =limitm2;
            end
        else
            m2begin=U(end,2);
            m2end  =limitm2;
        end
        [m2,ratio,peakm2,deltam]=onePeak(m2begin,m2end,firstRatio(pointPeakes(1)),EorO,U,fraction);
    elseif noumberOfPeakes>1
        ratio=[];
        m2=[];
        deltam=inf(1,noumberOfPeakes);
        peakm2=zeros(1,noumberOfPeakes);
        
        if pointPeakes(1)<pointsmin(1)
            for i=1:noumberOfPeakes
                if i==1
                    m2begin=U(end,2);
                    m2end  =first_m2Space(pointsmin(i));
                elseif i>numel(pointsmin)
                    m2begin=first_m2Space(pointsmin(i-1));
                    m2end=limitm2;
                else
                    m2begin=first_m2Space(pointsmin(i-1));
                    m2end  =first_m2Space(pointsmin(i));
                end
                fprintf('计算第%i个峰:\n', i);
                [m21,ratio1,peak_m2,deltam_i]=onePeak(m2begin,m2end,firstRatio(pointPeakes(i)),EorO,U,fraction);
                fprintf('m^2->：%15.10f\n', peak_m2);
                disp('------------------');
                m2=[m2 m21];
                ratio=[ratio ratio1];
                peakm2(i)=peak_m2;
                deltam(i)=deltam_i;
            end
        else
            for i=1:noumberOfPeakes
                if i==numel(pointsmin)
                    m2begin=first_m2Space(pointsmin(i));
                    m2end=limitm2;
                else
                    m2begin=first_m2Space(pointsmin(i));
                    m2end  =first_m2Space(pointsmin(i+1));
                end
                fprintf('计算第%i个峰:\n', i);
                [m21,ratio1,peak_m2,deltam_i]=onePeak(m2begin,m2end,firstRatio(pointPeakes(i)),EorO,U,fraction);
                fprintf('m^2->：%15.10f\n', peak_m2);
                disp('------------------');
                m2=[m2 m21];
                ratio=[ratio ratio1];
                peakm2(i)=peak_m2;
                deltam(i)=deltam_i;
            end
        end
        
        
       
    end 
    if m2end<limitm2  % add the last section.
        m2t =linspace(m2end,limitm2,200);
        ratiot=spectrum(EorO,U,m2t,200,fraction);
        m2=[m2 m2t];
        ratio=[ratio ratiot];
    end
    if  m2(1)>U(end,2)  % add the last section.
        
        m2t =linspace(U(end,2),m2(1),200);
        ratiot=spectrum(EorO,U,m2t,200,fraction);
        m2=[m2t m2];
        ratio=[ratiot ratio];
    end
    m=sqrt(m2);
    m2mDeltamTau=[peakm2' sqrt(peakm2)' deltam' 1./deltam'];
    
else
    % No peak found
    warning('No peak found!');
    m2=first_m2Space;
    ratio=firstRatio;
    m=sqrt(m2);
    m2mDeltamTau=[nan nan nan nan];

end

end%function



function [m2,ratio,peakm2,deltam]=onePeak(m2begin,limitm2,peakValue,EorO,U,fraction)
num=1000;%增加点数
peakValue1=peakValue;%保存原始峰值的大小
m2 =linspace(m2begin,limitm2,num);%这个区间内只有一个真正的峰。
ratio=spectrum(EorO,U,m2,num,fraction);%重新计算
[peakValue2,peakPostion]=max(ratio);%找出极大值点的位置，这次应该是一个峰值。
reltol=abs((peakValue2-peakValue1)/peakValue1);%确定相对误差
minm2=m2begin;
maxm2=limitm2;
i=1;%准备循环的标志的初值
%下面的while的循环式判断是否峰值满足要求。
[deltam2,fwhmYorN]=fwhm(m2,ratio);%deltam2 = delta m^2
[deltam]=fwhm2(m2,ratio);%deltam = delta m
while reltol>10-4%相对误差太大进如下一个循环
    [deltam2,fwhmYorN]=fwhm(m2,ratio);
	[deltam]=fwhm2(m2,ratio);
   % fprintf('fwhmYorN= %f\n',fwhmYorN);
    peakValue1=peakValue2;  
    if fwhmYorN==1              
        minm2t=minm2;% minm2t临时变量
        minm2=m2(peakPostion)-deltam2;%进一步缩小区间的范围
        if minm2<minm2t
            minm2=minm2t;
        end
        maxm2t=maxm2;% maxm2t 临时变量
        maxm2=m2(peakPostion) + deltam2;%
        if maxm2>maxm2t
            maxm2=maxm2t;
        end
      else
        num=num*2;
    end
    %
   % fprintf('num= %f\n',num)
    if num>=10^4 %防止死循环
        fprintf('循环结束\n')
        warning('Peak:if2','峰的性质没有满足要求可能需要增加循区间点的个数')
        fprintf('\n')
        break
    end
    m2 =linspace(minm2,maxm2,num);
    ratio=spectrum(EorO,U,m2,num,fraction);
    [peakValue2,peakPostion]=max(ratio);
    reltol=abs((peakValue2-peakValue1)/peakValue1); 
    fprintf('循环第 %d 次\n',i)   
    fprintf('峰值相对误差：reltol= %f\n',reltol)
    i=i+1;    
    if i>=10 %防止死循环
        fprintf('循环结束\n')
        warning('Peak:if2','峰的性质没有满足要求可能需要增加循环的次数')
        fprintf('\n')
        break
    end
    %      m2 =linspace(minm2t,maxm2t,num);
    %      ratio=spectrum(EorO,U,m2,num,fraction);
end
numfix=200;
if fwhmYorN==1
    fprintf('峰值为： %f\n',ratio(peakPostion));
    peakm2=m2(peakPostion);
	m2t1 =linspace(m2begin,min(m2),numfix);
	ratiot1=spectrum(EorO,U,m2t1,numfix,fraction); %求出没有共振峰的区间ratio的值
	m2t2 =linspace(max(m2),limitm2,numfix);
	ratiot2=spectrum(EorO,U,m2t2,numfix,fraction); %求出没有共振峰的区间ratio的值
	m2=[m2t1 m2 m2t2];%总m^2的区间加上没有共振峰的m^2的区间
    ratio=[ratiot1 ratio ratiot2];%加上没有共振峰的区间ratio的值
	else
    disp('这不是一个真正的峰，没有半高宽！');
	m2 =linspace(m2begin,limitm2,numfix);
    ratio=spectrum(EorO,U,m2,numfix,fraction);%重新计算;
	peakm2=nan;
    deltam2=inf;
	deltam=deltam2;
end
end%function   
    
%*****************************************************************************

%--------------------------------------------------------------------------
%给定m^2区间求出对应谱
function ratio=spectrum(EorO,U,m2,num,fraction)
%num：计算的步数,和m2矩阵的长度一致。
%m：几分之一
%是奇函数还是偶函数
if EorO==1
    %Odd:
    y10=0;
    y20=1;
else
    %Even:
    y10=1;
    y20=0;
end

x=U(:,1);
hx=x(2)-x(1);
U=U(:,2);
n=numel(U);
ratio = zeros(1,num);
for j=1:num
    V     = U-m2(j);
    y     = numerov(y10, y20, hx, V);
    dely1=sum(y.^2);
    dely2=sum(y(1,1:fix(n*fraction)).^2);%取1/10
    ratio(j)=dely2/dely1;
end

end%function
