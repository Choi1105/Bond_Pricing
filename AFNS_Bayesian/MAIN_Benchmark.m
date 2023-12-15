clear
clc

%% Options
lag = 2;                                       % 자료의 MA 시차 
tau = [3 6 9 12 18 24 30 36 60 120]';          % 만기 
maxac = 50;                                    % 비효율성 계수의 시차
freq = 100;                                    % 중간결과 
n1 = 25000;                                     % 샘플링 크기
n0 = 10000;                                     % 번인 크기

%% Data Loading
Macrom_0 = xlsread('KOR_YC_2021_08_M','Macro','B2:C249');
YCm_0 = xlsread('KOR_YC_2021_08_M','YC','B2:K249');

m = cols(Macrom_0);
% Moving Averaging order 3
% 계절성 or 노이즈로 인한 결과 왜곡을 방지하기 위해 smooth
% MA order가 너무 크면 너무 smooth해서 오히려 통계적으로 결과가 안좋은 경우 존재 
YCm = YCm_0;
Macrom = Macrom_0;
if lag > 0
    for t = lag+1:length(Macrom_0)
        mt = meanc(Macrom_0(t-lag:t, 1:m));
        Macrom(t, 1:m) = mt';
    end
end

% Macrom = detrend(Macrom);   % 추세제거시 사용
% 1200 나눠주는 이유 = 이론적인 모형과 부합
% Decimal and Monthly 
Macrom = 1*standdc(Macrom(lag+1:end, :))/1200;     % 표준화 or 평균제거
YCm = YCm(lag+1:end, :)/1200;

%% 추정
% 거시변수 설정 
Macrom_M = Macrom(:, 1); % FRP에 영향을 주는 Unspanned Risk 
Macrom_Z = Macrom(:, 2); % FRP에 영향을 주지 않는 Unspanned Risk

% 사전 설정
Spec = make_Spec(YCm, Macrom_M, Macrom_Z, tau, freq, maxac, n0, n1);

% 파라메터 샘플링
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
