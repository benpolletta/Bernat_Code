function Figure_4 % (hi_hr, cplot_norm)

hi_hr = 'drug'; cplot_norm = '';

phase_freqs = 1:.25:12; amp_freqs = 20:5:200;
no_pfs = length(phase_freqs); no_afs = length(amp_freqs);

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat')

stats={'median'}; %,'mean','std'}; 
no_stats = length(stats);
long_stats={'Median'}; % ,'Mean','St. Dev.'};

norms={''}; % , '_pct'}; 
no_norms = length(norms);
long_norms={''}; % , '% Change'};% From Baseline'};

no_pre=4; no_post=12;
[hr_labels, ~, long_hr_labels] = make_period_labels(no_pre, no_post, 'hrs');
no_hr_periods = length(hr_labels);

no_pre=4; no_post=16;
[BP_hr_labels, ~, ~] = make_period_labels(no_pre, no_post, 'hrs');
no_BP_hr_periods = length(BP_hr_labels);
short_BP_hr_labels = -4:16; short_BP_hr_labels(short_BP_hr_labels == 0) = [];

tick_spacing = floor(no_BP_hr_periods/5);

no_bands = 6;
            
bands_plotted = [2 5 6]; no_bands_plotted = length(bands_plotted);

drug_p_val_index = [1 4 2 3];

c_order = [0 0 1; 0 .5 0; 1 0 0];
    
clear titles xlabels ylabels
           
All_cplot_data = nan(no_afs, no_pfs, no_drugs, no_hr_periods, no_stats, no_channels, no_norms);
    
for n=1:no_norms
    
    for c=1:no_channels
        
        ch_dir = ['ALL_', channel_names{c}];
            
        %% Getting colorplot data.
        
        for d = 1:(no_drugs - 1)
            
            load([ch_dir,'/',ch_dir,'_p0.99_IEzs_MI','/',...
                ch_dir,'_p0.99_IEzs_hrMI_hr',norms{n},'_',drugs{d},'_cplot_data.mat'])
            
            All_cplot_data(:, :, d, :, :, c, n) = MI_stats(:, :, 1, :, 1);
            
        end
        
    end
    
end

if strcmp(hi_hr, 'independent')

    [~, max_hr_indices] = nanmax(nanmax(nanmax(All_cplot_data(:, :, :, 5:end, :, :, :))), [], 4);
    
    max_hr_indices = reshape(max_hr_indices, no_drugs, no_stats, no_channels, no_norms);

elseif strcmp(hi_hr, 'drug')
   
    [~, max_hr_indices] = nanmax(nanmax(nanmax(nanmax(All_cplot_data(:, :, :, 5:end, :, :, :))), [], 6), [], 4);
    
    max_hr_indices = repmat(reshape(max_hr_indices, no_drugs, no_stats, 1, no_norms), [1 1 no_channels 1]);
    
end

handle = nan(no_drugs, no_norms, no_stats);

All_cplot_for_plot = nan(no_afs, no_pfs, no_drugs + 1, no_stats, no_channels, no_norms);

for n = 1:no_norms
    
    for c = 1:no_channels
        
        for s = 1:no_stats
            
            for d = 1:(no_drugs - 1)
                
                All_cplot_for_plot(:, :, d, s, c, n) = All_cplot_data(:, :, d, max_hr_indices(d, s, c, n) + 4, s, c, n); %...
                    % - All_cplot_data(:, :, 1, max_hr_indices(d, s, c, n) + 4, s, c, n);
                
            end
            
        end
        
    end
    
end

if strcmp(cplot_norm, '_row')
    
    max_by_channel = reshape(nanmax(nanmax(nanmax(All_cplot_for_plot))), no_stats, no_channels, no_norms);
    
    min_by_channel = reshape(nanmin(nanmin(nanmin(All_cplot_for_plot))), no_stats, no_channels, no_norms);
    
elseif strcmp(cplot_norm, '_col')
    
    max_by_drug = reshape(nanmax(nanmax(nanmax(All_cplot_for_plot)), [], 5), no_drugs + 1, no_stats, no_norms);
    
    min_by_drug = reshape(nanmin(nanmin(nanmin(All_cplot_for_plot)), [], 5), no_drugs + 1, no_stats, no_norms);
    
end

n = 1; s = 1; figure;

for c = 1:no_channels
    
    for d = 2:(no_drugs - 1)
        
        %% Plotting comodulograms.
        
        subplot(no_channels, 2*(no_drugs - 2), (c - 1)*2*(no_drugs - 2) + d - 1)
        
        imagesc(phase_freqs, amp_freqs, All_cplot_for_plot(:, :, d, s, c, n))
        
        if c == 1 && d == 3
           
            hold on
            
            plot([1 4 4 1 1], [120 120 160 160 120], 'LineWidth', 3, 'Color', [0 1 0])
            
            plot([6 11 11 6 6], [120 120 160 160 120], 'LineWidth', 3, 'Color', [1 0 0])
            
        end
        
        if strcmp(cplot_norm, '_row')
            
            caxis([min_by_channel(s, c, n) max_by_channel(s, c, n)])
            
            if d == no_drugs - 1
                
                colorbar('FontSize', 16)
                
            end
            
        elseif strcmp(cplot_norm, '_col')
            
            caxis([min_by_drug(d, s, n) max_by_drug(d, s, n)])
            
            colorbar('FontSize', 16)
            
        else
            
            % c_lims = caxis;
            %
            % caxis([max(c_lims(1), -abs(c_lims(2))) c_lims(2)])
            
            colorbar('FontSize', 16)
            
        end
        
        axis xy
        
        set(gca, 'FontSize', 16)
        
        if c == 1
            
            title({drugs{d}; [long_hr_labels{max_hr_indices(d, s, c, n) + 4}]}) % long_stats{s}, ' MI, ', long_norms{n},
            
            % else
            %
            %     title(long_hr_labels{max_hr_indices(d - 1, s, c, n) + 4 - 1})
            
        end
        
        if d == 2, ylabel({channel_names{c}; 'Amp. Freq. (Hz)'}), end
        
        if c == no_channels, xlabel('Phase Freq. (Hz)'), end
        
    end
    
end

significance = .025;
sig_flag = make_label('p', significance, .025);
            
bands_plotted = 5:6;
no_bands_plotted = length(bands_plotted);
band_label = make_label('bands', bands_plotted, 1:6);

phase_freqs = 1:.25:12; amp_freqs = 20:5:200;
no_pfs = length(phase_freqs); no_afs = length(amp_freqs);

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat')

stats={'median'}; %,'mean','std'}; 
no_stats = length(stats);
long_stats={'Median'}; % ,'Mean','St. Dev.'};

norms={''}; % , '_pct'}; 
no_norms = length(norms);
long_norms={''}; % , '% Change'};% From Baseline'};

no_pre=4; no_post=12;
[BP_hr_labels, ~, ~] = make_period_labels(no_pre, no_post, 'hrs');
no_BP_hr_periods = length(BP_hr_labels);
short_BP_hr_labels = -no_pre:no_post; short_BP_hr_labels(short_BP_hr_labels == 0) = [];

no_pre=2; no_post=6;
[BP_6min_labels, ~, ~] = make_period_labels(no_pre, no_post, '6mins');
no_BP_6min_periods = length(BP_6min_labels);
short_BP_6min_labels = ((-no_pre*60 + 3):6:(no_post*60 - 3))/60; short_BP_6min_labels(short_BP_6min_labels == 0) = [];

timesteps = {'hr', '6min'};
timestep_labels = {short_BP_hr_labels, short_BP_6min_labels};
no_timesteps = length(timestep_labels);
timestep_lengths = cellfun(@(x) length(x), timestep_labels);
tick_spacing = floor(timestep_lengths/4);

no_bands = 6;

drug_p_val_index = [1 4 2 3];

c_order = distinguishable_colors(no_bands_plotted);
c_order(bands_plotted == 5, :) = [0 1 0];
c_order(bands_plotted == 6, :) = [1 0 0];

clear titles xlabels ylabels

All_BP_stats = nan(max(timestep_lengths), no_channels, no_bands, no_stats, no_drugs, no_norms);

load('summed_MI_ranksum_by_subject.mat')


n = 1; t = 2; s = 1;

suffix = ['_summed_hrMI', norms{n}, '_', timesteps{t}, '_BP_stats'];

for c=1:no_channels
    
    ch_dir = ['ALL_', channel_names{c}];
    
    %% Getting time series data.
    
    load([ch_dir, '/', ch_dir, suffix, '.mat'])
    
    BP_stats_new = permute(BP_stats, [2, 1, 3, 4]);
    
    BPs_dims = size(BP_stats_new);
    
    BP_stats_new = reshape(BP_stats_new, BPs_dims(1), 1, BPs_dims(2), BPs_dims(3), BPs_dims(4));
    
    All_BP_stats(1:timestep_lengths(t), c, :, :, :, n, t) = BP_stats_new(1:timestep_lengths(t), :, :, 1, :); % [1 4 5] are where the median, mean, and std are.
    
    % Getting band_labels.
    
    ranksum_suffix = ['hrMI', norms{n}, '_6min_ranksum'];
    
    load([ch_dir, '/', ch_dir, '_summed_', ranksum_suffix, '.mat'], 'band_labels')
    
end

% Bonferroni correcting & testing p-values.

All_BP_test(:, :, :, :, :, n, t) = All_BP_ranksum(:, :, :, :, :, n, t) <= significance;

for c = 1:no_channels
    
    %% Plotting time series.
    
    for d = 2:(no_drugs - 1)
        
        clear plot_stats plot_test
        
        plot_stats = squeeze(All_BP_stats(:, c, bands_plotted, s, d, n, t) - All_BP_stats(:, c, bands_plotted, s, 1, n, t));
        
        ax(c, d - 1) = subplot(no_channels, 2*(no_drugs - 2), (c - 1)*2*(no_drugs - 2) + 2 + d - 1);
        
        set(gca, 'NextPlot', 'add', 'ColorOrder', c_order)
        
        plot((1:timestep_lengths(t))', zeros(size(1:timestep_lengths(t)))', 'k--')
        
        hold on
        
        plot((1:timestep_lengths(t))', plot_stats(1:timestep_lengths(t), :), 'LineWidth', 2)
        
        axis tight
        
        box off
        
    end
    
    % sync_axes(ax(c, :))
    
    %% Plotting stats.
    
    for d = 2:(no_drugs - 1)
        
        subplot(no_channels, 2*(no_drugs - 2), (c - 1)*2*(no_drugs - 2) + 2 + d - 1)
        
        plot_test = double(All_BP_test(1:timestep_lengths(t), c, :, bands_plotted, d - 1, n, t)); %]; % drug_p_val_index(d) - 1)]; [nan(size(All_BP_test(:, :, bands_plotted(b), d - 1)))... % drug_p_val_index(d) - 1)))...
        
        plot_test = permute(squeeze(plot_test), [1 3 2]);
        
        plot_test = reshape(plot_test, size(plot_test, 1), no_bands_plotted*2);
        
        plot_test(plot_test == 0) = nan;
        
        side_vec = [ones(1, no_bands_plotted) zeros(1, no_bands_plotted)];
        
        add_stars(gca, (1:timestep_lengths(t))', plot_test(1:timestep_lengths(t), :), side_vec, c_order)
        
        set(gca, 'XTick', 1:tick_spacing(t):timestep_lengths(t), 'XTickLabel', round(timestep_labels{t}(1:tick_spacing(t):end)), 'FontSize', 16)
        
        if c == 1
            
            title(drugs{d})
            
            if d - 1 == 1
                
                legend(band_labels(bands_plotted), 'Location', 'NorthEast', 'FontSize', 6)
                
            end
            
        elseif c == no_channels
            
            xlabel('Time Rel. Inj. (h)')
            
        end
        
        if d - 1 == 1
            
            ylabel('MI (sum)')
            
        end
        
    end
    
    % sync_axes(ax(c, :))
    
end

save_as_eps(gcf, 'Figure_4') %, 'orientation', 'portrait')

% save_as_pdf(gcf, 'Figure_4')


end
