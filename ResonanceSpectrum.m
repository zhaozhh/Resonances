function ResonanceSpectrum
%2010��12��Ϊ����ComplexPhi4ConformalFlat��д��
%������������ӵ�����ƽ���ף��������ص����ݡ�
%flag:������ʶʱ�յ����ʣ���ֵ�ֱ�Ϊ1,0��-1��ӦdS, flat, AdSʱ�ա�
%Lambda4��4Dʱ�յ�����ѧ����ֵ��
%EorO:    ȷ�����̽����ż�ԣ�1���溯����2��ż������
%chiral�� ȷ�������ӵ������ԣ�1�����ַ����ӣ�2���ַ����ӡ�
%fraction:��С������������ı�ֵ��
%a:       һ���ɵ��Ĳ������ں���ComplexPhi4ConformalFlat�б�ʾ���Ǻ��¶��йص�
%         ����������������������
clc
pV=input('���ܺ������ļ�����Ĭ����xxx.dat����  ','s');
if isempty(pV)
    pV ='\URZ.dat';
end

pathname=pwd;
U=load([pathname,'\',pV]);

EorO=input('����������ƣ�1-�溯����2-ż���� (Ĭ��Ϊ1)��  ');
if isempty(EorO)
    EorO = 1;
end
fraction=input('������Ϊƽ�沨ʱ�ֶԼ��ʵ�ֵ��Ĭ��Ϊ0.1����  ');
if isempty(fraction)
    fraction =0.1;
end


%�������ݣ�
%�������ݣ�
if EorO==1
    str_EorO='_Odd';
else
    str_EorO='_Even';
end


%U=load(potentialName);
[m,m2,ratio,m2mDeltamTau]=findpeaks2(EorO,U,fraction);
%dlmwrite ([pathname,'\m2mDeltamTau',str_EorO,'.csv'],m2mDeltamTau, 'precision', '%.5f')
% save([pathname,'\m2mDeltamTau',str_EorO,'.dat'],'m2mDeltamTau','-ASCII','-double')
figure
plot(m2(1,2:numel(m2)),ratio(1,2:numel(m2)),'k-')
%    axis([0 max(U)*(1+0.2) 0 1])%�趨����ķ�Χ��
%    title(str3)
xlabel('m^2')
ylabel('P_R')
data=[m2(1,2:numel(m2));m(1,2:numel(m2));ratio(1,2:numel(m2))]';
% open a file for writing
path2=[pathname,'\m^2_m_Deltam_Tau',str_EorO,'.csv'];
m2mRelativeProbability=[pathname,'\m^2_m_RelativeProbability',str_EorO,'.dat'];
save(m2mRelativeProbability,'data','-ASCII','-double');
%���m^2,m,delta m,tau
fid = fopen(path2, 'w');
fprintf(fid, 'm^2, m, delta m , tau \n\n');
fprintf(fid, '%.5f,   %.5f, %.15f,  %.5f \n', m2mDeltamTau');
fclose(fid);
