function NVP_deltaHFO_thetaHFO_plot_by_state(drug, state, subject, hours)

if nargin < 3, subject = []; end
if isempty(subject), subject = 'A99'; end

if nargin < 4, hours = []; end
if isempty(hours), hours = 1:4; end
hours_label = sprintf('post%dto%d', min(hours), max(hours));

load('drugs.mat')

drug_indicator = strcmp(drugs, drug);

load('states.mat')

state_indicator = strcmp(states, state);

if nargin < 1, drug = ''; end
if isempty(drug), drug = 'NVP'; end

load('AP_freqs')

figure

%% Plotting Occipital MI for hour 4.
            
load(['ALL_Occipital/ALL_Occipital_p0.99_IEzs_MI/',...
    'ALL_Occipital_p0.99_IEzs_hrMI_4hr_by_state_cplot_data.mat'])
        
pd_indicator = strcmp(cat2_labels, hours_label);

subplot(2, 3, 1)

imagesc(phase_freqs, amp_freqs, MI_stats(:, :, state_indicator, pd_indicator, drug_indicator, 1))

axis xy

colorbar

title({'All Occipital PAC'; sprintf('Hrs. %d to %d Post-Injection', min(hours), max(hours));...
        [long_states{state_indicator}, ', Median MI']}, 'FontSize', 16)

ylabel('Amp. Freq. (Hz)', 'FontSize', 14)

xlabel('Phase Freq. (Hz)', 'FontSize', 14)

%% Plotting top quartiles for delta-HFO and theta-HFO PAC.

load(sprintf('dHFO_tHFO_quantiles_%s_q0.25_%s.mat', hours_label, state))

HFO_labels = {'\delta-HFO', '\theta-HFO'};

for lf = 1:2
    
    subplot(2, 3, 1 + lf)
    
    imagesc(phase_freqs, amp_freqs, reshape(median_qPAC(:, lf, drug_indicator, 2), no_afs, no_pfs))
    
    axis xy
    
    colorbar
    
    title({['Top Quartile for ', HFO_labels{lf}, ' PAC'];...
        sprintf('Hrs. %d to %d Post-Injection', min(hours), max(hours));...
        [long_states{state_indicator}, ', Median MI']}, 'FontSize', 16)
    
    xlabel('Phase Freq. (Hz)', 'FontSize', 14)
    
end

%% Plotting delta-HFO against theta-HFO.

ch_dir = 'ALL_Occipital'; measure = 'p0.99_IEzs'; q = .25;

sMI = load([ch_dir, '/', ch_dir, '_', measure, '_summed.mat']);
sMI = sMI.summed_MI;
drugs = text_read([ch_dir,'/',ch_dir, '_', measure,'_drugs.txt'],'%s');
subjects = text_read([ch_dir,'/',ch_dir, '_', measure,'_subjects.txt'],'%s');
hrs = text_read([ch_dir,'/',ch_dir, '_', measure,'_hrs.txt'],'%s');
sixmins = text_read([ch_dir,'/',ch_dir, '_', measure,'_6mins.txt'],'%d');
states = text_read([ch_dir,'/',ch_dir, '_', measure,'_states.txt'],'%s');

subj_index = strcmp(subjects, subject) & strcmp(drugs, drug) & strcmp(states, state);
            
for six = 1:40
    
    six_index = abs(sixmins - (six*6 - 3)) < eps & subj_index;
    
    subj_deltaHFO(six, 1) = nanmean(sMI(six_index, end));
    
    subj_thetaHFO(six, 1) = nanmean(sMI(six_index, end - 1));
    
end

subj_deltaHFO = nanzscore(subj_deltaHFO);

subj_deltaHFO_high = subj_deltaHFO;

subj_deltaHFO_high(subj_deltaHFO < quantile(subj_deltaHFO, 1 - q)) = nan;

subj_thetaHFO = nanzscore(subj_thetaHFO);

subj_thetaHFO_high = subj_thetaHFO;

subj_thetaHFO_high(subj_thetaHFO < quantile(subj_thetaHFO, 1 - q)) = nan;

h = subplot(2, 2, 3);

sixmin_labels = (1:40)*6 - 3;
    
set(h, 'NextPlot', 'add', 'ColorOrder', [0 1 0; 1 0 0])

plot(sixmin_labels', [subj_deltaHFO subj_thetaHFO], 'LineWidth', 2)

hold on

plot(sixmin_labels', [subj_deltaHFO_high subj_thetaHFO_high], 'o', 'LineWidth', 2)

legend({'\delta-HFO PAC','\theta-HFO PAC'})

plot(sixmin_labels', zeros(size(subj_deltaHFO)), '--k')

axis tight

box off

title({['Occipital \delta-HFO and \theta-HFO PAC Post-Injection, ', long_states{state_indicator}]; '(Representative Individual)'}, 'FontSize', 16)

ylabel('Mean MI (a.u.)', 'FontSize', 14)

xlabel('Time (h) Rel. Injection', 'FontSize', 14)

%% Plotting correlation of delta-HFO and theta-HFO PAC.

load(['Occipital_', drug, '_thetaHFO_deltaHFO_multidrug_', state, '.mat'])

subplot(2, 3, 6)

imagesc(c{1}, c{2}, n)

axis xy

xlabel('\delta-HFO PAC (a. u.)', 'FontSize', 14)

ylabel('\theta-HFO PAC (a. u.)', 'FontSize', 14)

hold on

plot(c{1}, p_fit(1)*c{1} + p_fit(2), 'w')

title({'Correlation of Post-Injection \delta-HFO and \theta-HFO PAC'; ['\rho = ', num2str(rho, '%.3g'), ', p = ', num2str(p, '%.3g')]}, 'FontSize', 16)

save_as_pdf(gcf, [drug, '_dHFO_tHFO_figure_', subject, '_', state])
