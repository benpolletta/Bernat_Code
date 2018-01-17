function BP_figure_stats_by_subject_plot_by_channel(significance, bands_plotted) % (hi_hr, cplot_norm)

if nargin < 1, significance = []; end
if isempty(significance), significance = .025; end
sig_flag = make_label('p', significance, .025);
            
if nargin < 2, bands_plotted = []; end
if isempty(bands_plotted), bands_plotted = 1:6; end % [1 2]; % [1 2 5];
no_bands_plotted = length(bands_plotted);
band_label = make_label('bands', bands_plotted, 1:6);

hi_hr = 'drug'; cplot_norm = '';

phase_freqs = 1:.25:12; amp_freqs = 20:5:200;
no_pfs = length(phase_freqs); no_afs = length(amp_freqs);

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat')

stats={'median'}; %,'mean','std'}; 
no_stats = length(stats);
long_stats={'Median'}; % ,'Mean','St. Dev.'};

norms={'', '_pct'}; 
no_norms = length(norms);
long_norms={'', '% Change'};% From Baseline'};

timesteps = {'_hrs', '_6min'};
no_timesteps = length(timesteps);

no_pre=4; no_post=12;
[hr_labels, ~, long_hr_labels] = make_period_labels(no_pre, no_post, 'hrs');
no_hr_periods = length(hr_labels);

no_pre=4; no_post=12;
[BP_hr_labels, ~, ~] = make_period_labels(no_pre, no_post, 'hrs');
no_BP_hr_periods = length(BP_hr_labels);
short_BP_hr_labels = -no_pre:no_post; short_BP_hr_labels(short_BP_hr_labels == 0) = [];

no_pre=2; no_post=6;
[BP_6min_labels, ~, ~] = make_period_labels(no_pre, no_post, '6mins');
no_BP_6min_periods = length(BP_6min_labels);
short_BP_6min_labels = ((-no_pre*60 + 3):6:(no_post*60 - 3))/60; short_BP_6min_labels(short_BP_6min_labels == 0) = [];

timestep_labels = {short_BP_hr_labels, short_BP_6min_labels};
timestep_lengths = cellfun(@(x) length(x), timestep_labels);
tick_spacing = floor(timestep_lengths/4);

no_bands = 6;

drug_p_val_index = [1 4 2 3];

no_skipped = 3;
c_order = distinguishable_colors(no_bands_plotted + no_skipped);
c_order = c_order((no_skipped + 1):end, :);
c_order(1, :) = [0 0 0]; c_order(4, :) = [0 1 1];
% c_order = [1 0 1; 0 1 1];

clear titles xlabels ylabels

preinj_data = nan(no_afs, no_pfs, no_stats, no_channels, no_norms);

All_BP_stats = nan(no_BP_6min_periods, no_channels, no_bands, no_stats, no_drugs, no_norms, no_timesteps);
    
load('BP_ranksum_by_subject')

for n=1:no_norms
    
    for t = 1:no_timesteps
        
        suffix = [norms{n}, timesteps{t}, '_BP_stats'];
        
        ranksum_suffix = [norms{n}, timesteps{t}, '_ranksum'];
        
        for c=1:no_channels
            
            ch_dir = ['ALL_', channel_names{c}];
            
            %% Getting time series data.
            
            load([ch_dir, '/', ch_dir, '_BP', suffix, '.mat'])
            
            BP_stats_new = permute(BP_stats, [2, 1, 3, 4]);
            
            BPs_dims = size(BP_stats_new);
            
            BP_stats_new = reshape(BP_stats_new, BPs_dims(1), 1, BPs_dims(2), BPs_dims(3), BPs_dims(4));
            
            temp_length = min(no_BP_6min_periods, BPs_dims(1));
            
            All_BP_stats(1:temp_length, c, :, :, :, n, t) = BP_stats_new(1:temp_length, :, :, 1, :); % [1 4 5] are where the median, mean, and std are.
            
            %% Getting band_labels.
            
            load([ch_dir, '/', ch_dir, '_BP', ranksum_suffix, '.mat'], 'band_labels')
            
        end
        
        % Bonferroni correcting & testing p-values.
        
        All_BP_test(:, :, :, :, :, n, t) = All_BP_ranksum(:, :, :, :, :, n, t) <= significance; % /(3*sum_all_dimensions(~isnan(All_BP_ranksum(:, :, :, :, :, n, t))));
        
    end
    
end

last_drug = no_drugs; % - 1;
first_drug = 1; % 2;
no_drugs_plotted = last_drug - first_drug + 1;

for n = 1:no_norms
    
    for t = 1:no_timesteps
        
        for s = 1:no_stats
            
            handle(n, s) = figure;
            
            %% Plotting time series.
            
            for c = 1:no_channels
                
                for d = first_drug:last_drug
                    
                    clear plot_stats
                    
                    if d == 1
                        
                        plot_stats = squeeze(All_BP_stats(:, c, bands_plotted, s, d, n, t));
                        
                    else
                    
                        plot_stats = squeeze(All_BP_stats(:, c, bands_plotted, s, d, n, t) - All_BP_stats(:, c, bands_plotted, s, 1, n, t));
                    
                    end
                        
                    ax(c, d - first_drug + 1) = subplot(no_channels, no_drugs_plotted, (c - 1)*no_drugs_plotted + d - first_drug + 1); % (d - 2)*no_bands_plotted + b)
                    
                    set(gca, 'NextPlot', 'add', 'ColorOrder', c_order) % {'--','-','*','*'}
                    
                    plot((1:timestep_lengths(t))', plot_stats(1:timestep_lengths(t), :), 'LineWidth', 2)
                    
                    hold on
                    
                    set(gca, 'XTick', 1:tick_spacing(t):timestep_lengths(t), 'XTickLabel', round(timestep_labels{t}(1:tick_spacing(t):end)), 'FontSize', 16)
                    
                    axis tight
                    
                    if c == 1
                        
                        title(drugs{d})
                        
                        if d - 1 == 1
                            
                            legend(band_labels(bands_plotted),...
                                'Location', 'NorthEast', 'FontSize', 6)
                            
                        end
                        
                    elseif c == no_channels
                        
                        xlabel('Time Rel. Inj. (h)')
                        
                    end
                    
                    if d - 1 == 1
                        
                        ylabel(channel_names{c})
                        
                    end
                    
                end
                
                % sync_axes(ax(c, :))
                
                %% Plotting stats.
                
                for d = 2:last_drug
                    
                    clear plot_test
                    
                    plot_test = squeeze(All_BP_test(:, c, :, bands_plotted, d - 1, n, t)); %]; % drug_p_val_index(d) - 1)]; [nan(size(All_BP_test(:, :, bands_plotted(b), d - 1)))... % drug_p_val_index(d) - 1)))...
                    
                    plot_test = permute(plot_test, [1 3 2]);
                    
                    plot_test = double(reshape(plot_test, size(plot_test, 1), no_bands_plotted*2));
                    
                    plot_test(plot_test == 0) = nan;
                    
                    side_vec = [ones(1, no_bands_plotted) zeros(1, no_bands_plotted)];
                    
                    subplot(no_channels, no_drugs_plotted, (c - 1)*no_drugs_plotted + d - first_drug + 1)
                    
                    add_stars(gca, (1:timestep_lengths(t))', plot_test(1:timestep_lengths(t), :), side_vec, c_order)
                    
                end
                
                % sync_axes(ax(c, :))
                
            end
            
            save_as_pdf(gcf,['ALL_BP', norms{n}, '_multichannel_',...
                timesteps{t}, make_label('bands', bands_plotted), sig_flag, '_stats_by_subject_plot_by_channel']) %, 'orientation', 'portrait')
            
        end
        
    end
    
end

end
