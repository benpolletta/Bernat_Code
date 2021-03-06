function summed_MI_p_vals_by_state_stats_by_subject_plot_by_channel(significance, bands_plotted)

if nargin < 1, significance = []; end
if isempty(significance), significance = .025; end
sig_label = make_label('p', significance, .025);
            
if nargin < 2, bands_plotted = []; end
if isempty(bands_plotted), bands_plotted = 1:6; end
no_bands_plotted = length(bands_plotted);
band_label = make_label('bands', bands_plotted, 1:6);

phase_freqs = 1:.25:12; amp_freqs = 20:5:200;
no_pfs = length(phase_freqs); no_afs = length(amp_freqs);

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat')

load('states.mat'), no_states = length(states);
long_states = {'Active Wake', 'Quiet Wake/non-REM Sleep', 'REM Sleep'};

stats={'median'}; %,'mean','std'}; 
no_stats = length(stats);
long_stats={'Median'}; % ,'Mean','St. Dev.'};

norms={'', '_pct'}; 
no_norms = length(norms);
long_norms={'', '%\Delta'};% From Baseline'};

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

c_order = distinguishable_colors(2*no_bands_plotted);
c_order = c_order(no_bands_plotted + (1:no_bands_plotted), :);
c_order(bands_plotted == 5, :) = [0 1 0];
c_order(bands_plotted == 6, :) = [1 0 0];

clear titles xlabels ylabels

All_BP_stats = nan(max(timestep_lengths), no_channels, no_bands, no_states, no_drugs, no_norms, no_timesteps);

load('summed_MI_ranksum_by_subject_by_state.mat')

for n=1:no_norms
    
    for t = 1:no_timesteps
        
        suffix = ['_summed_hrMI', norms{n}, '_', timesteps{t}, '_by_state'];
        
        for c=1:no_channels
            
            ch_dir = ['ALL_', channel_names{c}];
            
            %% Getting time series data.
            
            load([ch_dir, '/', ch_dir, suffix, '.mat'])
            
            BP_stats_new = permute(BP_stats, [3, 1, 2, 5, 4]);
            
            BPs_dims = size(BP_stats_new);
            
            BP_stats_new = reshape(BP_stats_new, [BPs_dims(1), 1, BPs_dims(2:end)]);
            
            All_BP_stats(1:timestep_lengths(t), c, :, :, :, n, t) = BP_stats_new(1:timestep_lengths(t), :, :, :, :, 1); % [1 4 5] are where the median, mean, and std are.
            
            % Getting band_labels.
            
            ranksum_suffix = ['hrMI', norms{n}, '_6min_ranksum'];
            
            load([ch_dir, '/', ch_dir, '_summed_', ranksum_suffix, '.mat'], 'band_labels')
            
        end
        
    end
    
end

% Making legend.

bands_plotted_labels = band_labels(bands_plotted);
my_legend1 = cellfun(@(x) strcat('saline ', x), bands_plotted_labels, 'UniformOutput', false);
my_legend2 = cellfun(@(x) strcat('drug ', x), bands_plotted_labels, 'UniformOutput', false);
my_legend = {my_legend1{:}, my_legend2{:}};
        
% Testing p-values.

All_BP_test = All_BP_ranksum; % <= .025;

for n = 1:no_norms
    
    for t = 1:no_timesteps
        
        for state = 1:no_states
            
            for s = 1:no_stats
                
                handle(n, s, t) = figure;
                
                for c = 1:no_channels
                    
                    %% Plotting time series.
                    
                    for d = 2:no_drugs
                        
                        clear plot_stats plot_test
                        
                        plot_stats = squeeze(All_BP_stats(:, c, bands_plotted, state, 1, n, t));
                        % [squeeze(All_BP_stats(:, c, bands_plotted, s, 1, n, t)) squeeze(All_BP_stats(:, c, bands_plotted, s, d, n, t))];
                        % squeeze(All_BP_stats(:, c, bands_plotted, s, d, n, t) - All_BP_stats(:, c, bands_plotted, s, 1, n, t));
                        
                        ax(c, d - 1) = subplot(no_channels, no_drugs - 1, (c - 1)*(no_drugs - 1) + d - 1);
                        
                        set(gca, 'NextPlot', 'add', 'ColorOrder', c_order) % , 'LineStyleOrder', {':', '-'})
                        
                        plot((1:timestep_lengths(t))', plot_stats(1:timestep_lengths(t), :), 'LineWidth', 1)
                        
                        hold on
                        
                        plot_stats = squeeze(All_BP_stats(:, c, bands_plotted, state, d, n, t));
                        
                        plot((1:timestep_lengths(t))', plot_stats(1:timestep_lengths(t), :), 'LineWidth', 2.25)
                        
                        axis tight
                        
                        box off
                        
                    end
                    
                    % sync_axes(ax(c, :))
                    
                    %% Plotting stats.
                    
                    for d = 2:no_drugs
                        
                        subplot(no_channels, no_drugs - 1, (c - 1)*(no_drugs - 1) + d - 1)
                        
                        plot_test = double(All_BP_test(1:timestep_lengths(t), c, 1, bands_plotted, d - 1, n, t, state)); %]; % drug_p_val_index(d) - 1)]; [nan(size(All_BP_test(:, :, bands_plotted(b), d - 1)))... % drug_p_val_index(d) - 1)))...
                        
                        plot_test = squeeze(plot_test); % permute(squeeze(plot_test), [1 3 2]);
                        
                        % plot_test = reshape(plot_test, size(plot_test, 1), no_bands_plotted*2);
                        
                        plot_test(plot_test == 0) = nan;
                        
                        hold off
                        
                        plot((1:timestep_lengths(t))', plot_test)
                        
                        % side_vec = [ones(1, no_bands_plotted) zeros(1, no_bands_plotted)];
                        
                        % add_stars(gca, (1:timestep_lengths(t))', plot_test, side_vec, c_order)
                        
                        set(gca, 'XTick', 1:tick_spacing(t):timestep_lengths(t), 'XTickLabel', round(timestep_labels{t}(1:tick_spacing(t):end)), 'FontSize', 16)
                        
                        if c == 1
                            
                            title([long_states{state}, ', ', drugs{d}])
                            
                            if d - 1 == 1
                                
                                legend(my_legend2, 'Location', 'NorthEast', 'FontSize', 6)
                                
                            end
                            
                        elseif c == no_bands_plotted
                            
                            xlabel('Time Rel. Inj. (h)')
                            
                        end
                        
                        if d - 1 == 1
                            
                            ylabel(channel_names{c})
                            
                        end
                        
                    end
                    
                    % sync_axes(ax(c, :))
                    
                end
                
                save_as_pdf(gcf,['ALL_summed_MI',norms{n},'_multichannel_', timesteps{t}, band_label, sig_label, '_', states{state}, '_stats_by_subject_plot_by_channel_p_vals']) %, 'orientation', 'portrait')
                
            end
            
        end
        
    end
    
end
