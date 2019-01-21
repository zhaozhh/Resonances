function [m,m2,ratio,m2mDeltamTau]=findpeaks5(EorO,U,fraction)
%�����еķ�ֵ�����յ��׺���
%2014.01.19,�Ż���һ�£��ṩ��fwhm��������һ�������ķ塣
%2014.01.16,�����������ʼֵ��Ϊ���ܵı߽�ֵ��
%2013.08.13,�޸ģ�������û�ҵ���ᱨ��ֹͣ���е����⡣
%2013.03.19,�Ľ�������ÿ����ʵ�ķ嶼�а�߿�����ԡ�
%2010��12��Ϊ����ferimionSpectrum,Zhen-Hua Zhao, zhaozhh02@gmail.com��
%EorO:    ȷ�����̽����ż�ԣ�1���溯����2��ż������
%chiral�� ȷ�������ӵ������ԣ�1�����ַ����ӣ�2���ַ����ӡ�
%fraction:��С������������ı�ֵ��
%mm2Deltam:������Ӧ���λ�õ�m m^2�ͷ�İ�߿�deltam�ľ���
%ratio:С����Ĳ�����ƽ���Ļ��ֳ��Դ�����Ĳ�����ƽ���Ļ���
%m��   ��������ֵ
%m2��  ������ƽ����Ҳ����Ѧ���̷��̵ı���ֵ��
%npeaks=5;%��ĸ����Ĺ���ֵ��
%******************************************************
%��ȫm^2���������ҳ���ĸ���
%******************************************************
global addition


if EorO==1   
    fprintf('\n �溯������\n ');
elseif EorO==2
    fprintf('\n ż��������\n ');
else
    error(' Wrong with EorO\n ')    
end
fprintf('------------------\n');
fprintf('The maximum of the potential: %f \n', max(U(:,2)));
addition=max(U(:,2))*0.5;%������m^2�����ֵҪ���Ǻ��������ֵ��addition��
fprintf('\n m^2�����ֵҪ�����ܵ����ֵ��%.3f   \n ',addition);
fprintf('m�����ֵ������   ��%f \n', sqrt(max(U(:,2))+addition));
fprintf(' \n');
num=input('�����ʼɨ��ļ����(Ĭ��Ϊ500):   ');
if isempty(num)
    num =500;
end

fprintf('��ʼɨ��ļ����Ϊ :%i \n', num)
limitm2=max(U(:,2))+addition;
first_m2Space =linspace(U(end,2),limitm2,num);
firstRatio=spectrum(EorO,U,first_m2Space,num,fraction);
pointPeakes=find(diff(sign(diff(firstRatio)))<0)+1;%find out all the maximum.
fprintf('���ܵķ�ĸ���Ϊ  ��%i\n', length(pointPeakes));
disp('��Լλ���ڣ�');
fprintf('m^2->��%f\n', first_m2Space(pointPeakes));
disp('------------------');

pointsmin=find(diff(sign(diff(firstRatio)))>0)+1;%find out all the minimum.

%fprintf('pointsmin    ��%f \n', numel(pointsmin));

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
                fprintf('�����%i����:\n', i);
                [m21,ratio1,peak_m2,deltam_i]=onePeak(m2begin,m2end,firstRatio(pointPeakes(i)),EorO,U,fraction);
                fprintf('m^2->��%15.10f\n', peak_m2);
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
                fprintf('�����%i����:\n', i);
                [m21,ratio1,peak_m2,deltam_i]=onePeak(m2begin,m2end,firstRatio(pointPeakes(i)),EorO,U,fraction);
                fprintf('m^2->��%15.10f\n', peak_m2);
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
num=1000;%���ӵ���
peakValue1=peakValue;%����ԭʼ��ֵ�Ĵ�С
m2 =linspace(m2begin,limitm2,num);%���������ֻ��һ�������ķ塣
ratio=spectrum(EorO,U,m2,num,fraction);%���¼���
[peakValue2,peakPostion]=max(ratio);%�ҳ�����ֵ���λ�ã����Ӧ����һ����ֵ��
reltol=abs((peakValue2-peakValue1)/peakValue1);%ȷ��������
minm2=m2begin;
maxm2=limitm2;
i=1;%׼��ѭ���ı�־�ĳ�ֵ
%�����while��ѭ��ʽ�ж��Ƿ��ֵ����Ҫ��
[deltam2,fwhmYorN]=fwhm(m2,ratio);%deltam2 = delta m^2
[deltam]=fwhm2(m2,ratio);%deltam = delta m
while reltol>10-4%������̫�������һ��ѭ��
    [deltam2,fwhmYorN]=fwhm(m2,ratio);
	[deltam]=fwhm2(m2,ratio);
   % fprintf('fwhmYorN= %f\n',fwhmYorN);
    peakValue1=peakValue2;  
    if fwhmYorN==1              
        minm2t=minm2;% minm2t��ʱ����
        minm2=m2(peakPostion)-deltam2;%��һ����С����ķ�Χ
        if minm2<minm2t
            minm2=minm2t;
        end
        maxm2t=maxm2;% maxm2t ��ʱ����
        maxm2=m2(peakPostion) + deltam2;%
        if maxm2>maxm2t
            maxm2=maxm2t;
        end
      else
        num=num*2;
    end
    %
   % fprintf('num= %f\n',num)
    if num>=10^4 %��ֹ��ѭ��
        fprintf('ѭ������\n')
        warning('Peak:if2','�������û������Ҫ�������Ҫ����ѭ�����ĸ���')
        fprintf('\n')
        break
    end
    m2 =linspace(minm2,maxm2,num);
    ratio=spectrum(EorO,U,m2,num,fraction);
    [peakValue2,peakPostion]=max(ratio);
    reltol=abs((peakValue2-peakValue1)/peakValue1); 
    fprintf('ѭ���� %d ��\n',i)   
    fprintf('��ֵ�����reltol= %f\n',reltol)
    i=i+1;    
    if i>=10 %��ֹ��ѭ��
        fprintf('ѭ������\n')
        warning('Peak:if2','�������û������Ҫ�������Ҫ����ѭ���Ĵ���')
        fprintf('\n')
        break
    end
    %      m2 =linspace(minm2t,maxm2t,num);
    %      ratio=spectrum(EorO,U,m2,num,fraction);
end
numfix=200;
if fwhmYorN==1
    fprintf('��ֵΪ�� %f\n',ratio(peakPostion));
    peakm2=m2(peakPostion);
	m2t1 =linspace(m2begin,min(m2),numfix);
	ratiot1=spectrum(EorO,U,m2t1,numfix,fraction); %���û�й���������ratio��ֵ
	m2t2 =linspace(max(m2),limitm2,numfix);
	ratiot2=spectrum(EorO,U,m2t2,numfix,fraction); %���û�й���������ratio��ֵ
	m2=[m2t1 m2 m2t2];%��m^2���������û�й�����m^2������
    ratio=[ratiot1 ratio ratiot2];%����û�й���������ratio��ֵ
	else
    disp('�ⲻ��һ�������ķ壬û�а�߿�');
	m2 =linspace(m2begin,limitm2,numfix);
    ratio=spectrum(EorO,U,m2,numfix,fraction);%���¼���;
	peakm2=nan;
    deltam2=inf;
	deltam=deltam2;
end
end%function   
    
%*****************************************************************************

%--------------------------------------------------------------------------
%����m^2���������Ӧ��
function ratio=spectrum(EorO,U,m2,num,fraction)
%num������Ĳ���,��m2����ĳ���һ�¡�
%m������֮һ
%���溯������ż����
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
    dely2=sum(y(1,1:fix(n*fraction)).^2);%ȡ1/10
    ratio(j)=dely2/dely1;
end

end%function
