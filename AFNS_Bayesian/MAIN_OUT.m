clear
clc

load post_param
load Spec
load Post_FM
load Risk_time_varying_Macro
load Risk_time_varying_Latent
load Risk_tsm_3Dm
load Risk_tsm_3D
load Risk_const
load Risk_time_varying

%% Latent Factors
F1 = 1200*meanc(Post_FM(:,:,1));
F2 = 1200*meanc(Post_FM(:,:,2));
F3 = 1200*meanc(Post_FM(:,:,3));

YCm = Spec.YCm;
tau = Spec.tau;
ntau = rows(tau);

Proxym(:, 1) = YCm(:, end)*1200;
Proxym(:, 2) = (YCm(:, 1) - YCm(:, end))*1200;
Proxym(:, 3) = (2*YCm(:, 5) - YCm(:, 1) - YCm(:, end))*1200;

n = rows(YCm);
datat = 1:n;
xtick = 23:24:n;
xticklabel = {'Jan 03', 'Jan 05','Jan 07','Jan 09','Jan 11', 'Jan 13','Jan 15','Jan 17','Jan 19','Jan 21'};

y = F1;
proxy = demeanc(Proxym(:, 1));
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2.6 scrsz(4)/1.2])
subplot(3,1,1)
h = plot(datat, y, 'b-',datat, proxy, 'k--');
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc(y) + 1;
lowy = minc(y) - 1;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
legend('Level', 'Long rate','Location','Northeast')
title('(a) Level')

y = F2;
proxy = demeanc(Proxym(:, 2));
subplot(3,1,2)
h = plot(datat, y, 'b-',datat, proxy, 'k--');
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc(y) + 1;
lowy = minc(y) - 1;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
legend('Slope', 'Spread','Location','Southeast')
title('(b) Slope')

y = F3;
proxy = demeanc(Proxym(:, 3));
subplot(3,1,3)
h = plot(datat, y, 'b-',datat, proxy, 'k--');
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc(y) + 1;
lowy = minc(y) - 1;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
legend('Curvature', 'Difference between Spreads','Location','Northeast')
title('(c) Curvature')

%% Term Premium
Z = Risk_tsm_3D;

yticklabel = tau;
y = tau;
ytick = [6, 36, 60, 120]';

uppz = maxc1(Z);
n =  rows(Z);
ztick = 0:1:uppz;

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/1.2 scrsz(3)/1.2 scrsz(4)/1.2])
subplot(2,1,1)
surf(datat, y, Z')
xlim([0 n+1])
zlim([-0.5 1.1*uppz])
xlabel('Time')
ylabel('Maturity')
zlabel('(%)')
set(gca,'ZTick', ztick);
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel);
set(gca,'XTickLabel', xticklabel,'FontSize', 9);
set(gca,'YTick', ytick);
set(gca,'YTickLabel');
title('Term Structure of Term Premium (3D)')

subplot(2,1,2)
y = Risk_tsm_3D;
h = plot(datat, y);
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time')
ylabel('(%)')
uppy = maxc1(y) + 0.5;
lowy = minc1(y) - 0.5;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
title('Term Structure of Term Premium (2D)')

%% Term Premium Decomposition
RP = Risk_tsm_3D;
y = RP;
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2.5 scrsz(3)/2.5 scrsz(4)/1.2])
subplot(4, 1, 1)
h = plot(datat, y);
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc1(y) + 0.5;
lowy = minc1(y) - 0.5;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
title('(a) Term premiums') 

y = Risk_const;
subplot(4, 1, 2)
h = plot(datat, y);
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc1(y) + 0.5;
lowy = minc1(y) - 0.5;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
title('(a) Time-invariant component') 

y = Risk_time_varying_Latent;
subplot(4, 1, 3)
h = plot(datat, y);
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc1(y) + 0.5;
lowy = minc1(y) - 0.5;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
title('(c) Time-varying component (Latent)') 

y = Risk_time_varying_Macro;
subplot(4, 1, 4)
h = plot(datat, y);
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
xlabel('Time','FontSize',9)
ylabel('(%)','FontSize',9)
uppy = maxc1(y) + 0.5;
lowy = minc1(y) - 0.5;
ylim([lowy uppy])
xlim([0 n+1])
set(h,'LineWidth', 1.5)
title('(d) Time-varying component (Macro)') 

%% Term Premium and EH with Credible Interval
tau_j = 10; 
spread = 1200*(YCm(:,tau_j) - YCm(:,1));
Risk_tsm_3Dm_taujm = Risk_tsm_3Dm(:, :, tau_j);
p = [0.05 0.5 0.95];
Risk_tsm_3Dm_tauj = quantile(Risk_tsm_3Dm_taujm,p)';
EHm_tauj = kron(ones(rows(Risk_tsm_3Dm_taujm),1), spread') - Risk_tsm_3Dm_taujm;
EH_tauj = quantile(EHm_tauj,p)';

figure
subplot(2,1,1)
h=plot(datat, Risk_tsm_3Dm_tauj(:, 1), 'b--',datat, Risk_tsm_3Dm_tauj(:, 2), 'k-',datat, Risk_tsm_3Dm_tauj(:, 3), 'b--' );
legend('5%', 'mean', '95%')
ylabel('(%)')
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
set(h,'LineWidth', 1.5)
title('Term Premium (10-year)')

subplot(2,1,2)
h=plot(datat, EH_tauj(:, 1), 'b--',datat, EH_tauj(:, 2), 'k-',datat, EH_tauj(:, 3), 'b--', datat, zeros(n,1), 'k:');
legend('5%', 'mean', '95%')
ylabel('(%)')
set(gca,'XTick', xtick);
set(gca,'XTickLabel', xticklabel);
set(gca,'XTickLabel', xticklabel,'FontSize',9);
set(h,'LineWidth', 1.5)
title('Expectation Hypothesis (10-year) - Short Rates')
