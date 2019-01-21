function WaveFunction
%EorO=:���溯������ż������1���溯����2��������
EorO=input('����������ƣ�1-�溯����2-ż���� (Ĭ��Ϊ1)��  ');
if isempty(EorO)
    EorO = 1;
end
Cm2=input('m^2��  ');
if isempty(Cm2)
    Cm2 =105.900363150622;
end
pV=input('���ܺ������ļ�����Ĭ����xxx.dat����  ','s');
if isempty(pV)
    pV ='\xxx.dat';
end

%�������ݣ�
%�������ݣ�


pathname=pwd;
U=importdata([pathname,'\',pV]);




m2  = Cm2;

if EorO==1
%Odd:
    y10=0;
    y20=1;
elseif EorO==2
%Even:
    y10=1;
    y20=0;
else
    fprint('****************************\n')
    fprint('The value of EorO is Wrong.\n')
    fprint('****************************\n')
end
%-----------


    U=U';
    x=U(1,:);
    hm2=x(2)-x(1);
    V=U(2,:);
    %������ܵ�ͼ��
    figure
    plot(x,V)
    xlabel('z')
    ylabel('U_{Lf}')
    %�󲨺���
    fprintf('Begin.\n');
    tic;%��������ʱ�䣬��ʼ
      V     = V-m2;
      y     = numerov(y10, y20, hm2, V);
    fprintf('The height of potential U is  %f \n', max(U(2,:)))
    %�������
    n=numel(x);
    x0=x;
    y1=y(1,:);
    x0=x0(n:-1:1);
    y1=y1(n:-1:1);
    x = [x0*(-1) x ];
    if EorO==1
    %Odd:
        y = [y1*(-1) y(1,:)];
        data=[x;y]';
        save([pathname,'\OddFunction.dat'],'data','-ASCII')       
    elseif EorO==2
    %Even:
        y = [y1 y(1,:)];
        data=[x;y]';
        save([pathname,'\EvenFuncion.dat'],'data','-ASCII')
    else
        fprint('****************************\n')
        fprint('The value of EorO is Wrong.\n')
        fprint('****************************\n')
    end    
    %�����������ͼ��
    figure%
    plot(x,y,'k-')
    if Cm2==0 %�������ģ
        axis([-10 10 -1 1])%�趨����ķ�Χ��
    end
    m2_str=num2str(m2);
    m2_str=['m^2= ' m2_str];
    title(m2_str)
    %axis([-10 10 0 1])
    xlabel('z')
    ylabel('\psi_f')
    data=[x; y]';
% open a file for writing

wave=[pathname,'\wavefunction_',m2_str,'.dat'];
save(wave,'data','-ASCII','-double');
