clear
clc

%% Options
lag = 2;                                       % �ڷ��� MA ���� 
tau = [3 6 9 12 18 24 30 36 60 120]';          % ���� 
maxac = 50;                                    % ��ȿ���� ����� ����
freq = 100;                                    % �߰���� 
n1 = 25000;                                     % ���ø� ũ��
n0 = 10000;                                     % ���� ũ��

%% Data Loading
Macrom_0 = xlsread('KOR_YC_2021_08_M','Macro','B2:C249');
YCm_0 = xlsread('KOR_YC_2021_08_M','YC','B2:K249');

m = cols(Macrom_0);
% Moving Averaging order 3
% ������ or ������� ���� ��� �ְ��� �����ϱ� ���� smooth
% MA order�� �ʹ� ũ�� �ʹ� smooth�ؼ� ������ ��������� ����� ������ ��� ���� 
YCm = YCm_0;
Macrom = Macrom_0;
if lag > 0
    for t = lag+1:length(Macrom_0)
        mt = meanc(Macrom_0(t-lag:t, 1:m));
        Macrom(t, 1:m) = mt';
    end
end

% Macrom = detrend(Macrom);   % �߼����Ž� ���
% 1200 �����ִ� ���� = �̷����� ������ ����
% Decimal and Monthly 
Macrom = 1*standdc(Macrom(lag+1:end, :))/1200;     % ǥ��ȭ or �������
YCm = YCm(lag+1:end, :)/1200;

%% ����
% �Žú��� ���� 
Macrom_M = Macrom(:, 1); % FRP�� ������ �ִ� Unspanned Risk 
Macrom_Z = Macrom(:, 2); % FRP�� ������ ���� �ʴ� Unspanned Risk

% ���� ����
Spec = make_Spec(YCm, Macrom_M, Macrom_Z, tau, freq, maxac, n0, n1);

% �Ķ���� ���ø�
[theta_st, post_param, post_summary, lnpostm, lnlikm, Post_FM, Spec] = MCMC_AFNS_SML(Spec);

save post_param post_param
save post_summary post_summary
save lnpostm lnpostm
save lnlikm lnlikm
save Spec Spec
save Post_FM Post_FM

[Risk_tsm_3Dm, Risk_tsm_3D, Risk_const, Risk_time_varying, Risk_time_varying_Latent, Risk_time_varying_Macro] = ...
    Term_premium(post_param, Spec);

save Risk_time_varying_Macro Risk_time_varying_Macro
save Risk_time_varying_Latent Risk_time_varying_Latent
save Risk_tsm_3Dm Risk_tsm_3Dm
save Risk_tsm_3D Risk_tsm_3D
save Risk_const Risk_const
save Risk_time_varying Risk_time_varying
