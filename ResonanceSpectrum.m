function ResonanceSpectrum
%2010年12月为函数ComplexPhi4ConformalFlat而写。
%功能是求费米子的质量平方谱，即输出相关的数据。
%flag:用来标识时空的性质，其值分别为1,0，-1对应dS, flat, AdS时空。
%Lambda4：4D时空的宇宙学常数值。
%EorO:    确定方程解得奇偶性，1：奇函数，2：偶函数。
%chiral： 确定费米子的手征性，1，左手费米子，2右手费米子。
%fraction:最小区间和最大区间的比值。
%a:       一个可调的参数，在函数ComplexPhi4ConformalFlat中表示的是和温度有关的
%         背景标量场的质量参数。
clc
%pV=input('势能函数的文件名（默认是xxx.dat）：  ','s');
%if isempty(pV)
    pV ='\URZ.dat';
%end

pathname=pwd;
U=load([pathname,'\',pV]);

%EorO=input('波函数的宇称，1-奇函数，2-偶函数 (默认为1)：  ');
%if isempty(EorO)
    EorO = 1;
%end
%fraction=input('波函数为平面波时现对几率的值（默认为0.1）：  ');
%if isempty(fraction)
    fraction =0.1;
%end


%导入数据：
%导入数据：
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
%    axis([0 max(U)*(1+0.2) 0 1])%设定坐标的范围。
%    title(str3)
xlabel('m^2')
ylabel('P_R')
data=[m2(1,2:numel(m2));m(1,2:numel(m2));ratio(1,2:numel(m2))]';
% open a file for writing
path2=[pathname,'\m^2_m_Deltam_Tau',str_EorO,'.csv'];
m2mRelativeProbability=[pathname,'\m^2_m_RelativeProbability',str_EorO,'.dat'];
save(m2mRelativeProbability,'data','-ASCII','-double');
%输出m^2,m,delta m,tau
fid = fopen(path2, 'w');
fprintf(fid, 'm^2, m, delta m , tau \n\n');
fprintf(fid, '%.5f,   %.5f, %.15f,  %.5f \n', m2mDeltamTau');
fclose(fid);
