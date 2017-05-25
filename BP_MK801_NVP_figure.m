function BP_MK801_NVP_figure % (hi_hr, cplot_norm)

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
            
bands_plotted = [1 2 5]; no_bands_plotted = length(bands_plotted);

drug_p_val_index = [1 4 2 3];

c_order = [0 0 1; 0 .5 0; 1 0 0];
    
clear titles xlabels ylabels

preinj_data = nan(no_afs, no_pfs, no_stats, no_channels, no_norms);

All_BP_stats = nan(no_BP_6min_periods, no_channels, no_bands, no_stats, no_drugs, no_norms, no_timesteps);

[All_BP_ranksum, All_BP_test] = deal(nan(no_BP_6min_periods, no_channels, 2, no_bands, no_drugs - 1, no_norms, no_timesteps));
    
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
            
            %% Getting stats p-values.
            
            load([ch_dir, '/', ch_dir, '_BP', ranksum_suffix, '.mat'])
            
            BP_ranksum_new = permute(BP_ranksum, [2, 4, 1, 3]);
            
            BPr_dims = size(BP_ranksum_new);
            
            temp_length = min(no_BP_6min_periods, BPr_dims(1));
            
            BP_ranksum_new = reshape(BP_ranksum_new, BPr_dims(1), 1, BPr_dims(2), BPr_dims(3), BPr_dims(4));
            
            All_BP_ranksum(1:temp_length, c, :, :, :, n, t) = BP_ranksum_new(1:temp_length, :, :, :, :);
            
        end
        
        % Bonferroni correcting & testing p-values.
        
        All_BP_test(:, :, :, :, :, n, t) = All_BP_ranksum(:, :, :, :, :, n, t) <= .01/(3*sum_all_dimensions(~isnan(All_BP_ranksum(:, :, :, :, :, n, t))));
        
    end
    
end

for n = 1:no_norms
    
    for t = 1:no_timesteps
        
        for s = 1:no_stats
            
            handle(n, s) = figure;
            
            %% Plotting time series w/ stats.
            
            for b = 1:no_bands_plotted
                
                for d = 2:(no_drugs - 1)
                    
                    clear plot_stats
                    
                    plot_stats = All_BP_stats(:, :, bands_plotted(b), s, d, n, t) - All_BP_stats(:, :, bands_plotted(b), s, 1, n, t);
                    % [All_BP_stats(:, :, bands_plotted(b), s, 1) All_BP_stats(:, :, bands_plotted(b), s, d)];
                    
                    ax(b, d - 1) = subplot(no_bands_plotted, no_drugs - 2, (b - 1)*(no_drugs - 2) + d - 1); % (d - 2)*no_bands_plotted + b)
                    
                    set(gca, 'NextPlot', 'add', 'LineStyleOrder', {'-','*','*'}, 'ColorOrder', c_order) % {'--','-','*','*'}
                    
                    plot((1:timestep_lengths(t))', plot_stats(1:timestep_lengths(t), :))
                    
                    hold on
                    
                    set(gca, 'XTick', 1:tick_spacing(t):timestep_lengths(t), 'XTickLabel', timestep_labels{t}(1:tick_spacing(t):end), 'FontSize', 16)
                    
                    axis tight
                    
                    if b == 1
                        
                        title(drugs{d})
                        
                        if d - 1 == 1
                            
                            legend({'Fr., drug - sal.', 'Occi., drug - sal.', 'CA1, drug - sal.'},...
                                'Location', 'NorthEast', 'FontSize', 6)
                            
                        end
                        
                    elseif b == no_bands_plotted
                        
                        xlabel('Time Rel. Inj. (h)')
                        
                    end
                    
                    if d - 1 == 1
                        
                        ylabel(band_labels{bands_plotted(b)})
                        
                    end
                    
                end
                
                linkaxes(ax(b, :))
                
                if b == no_bands_plotted
                    
                    y_lims = ylim;
                    
                    y_lim_upper = y_lims(2);
                    
                    axis(ax(b, 2), 'tight')
                    
                    y_lims = ylim(ax(b, 2));
                    
                    ylim([y_lims(1) y_lim_upper])
                    
                end
                
                linkaxes(ax(b, :), 'off')
                
                for d = 2:(no_drugs - 1)
                    
                    clear plot_test
                    
                    plot_test = All_BP_test(:, :, :, bands_plotted(b), d - 1, n, t); %]; % drug_p_val_index(d) - 1)]; [nan(size(All_BP_test(:, :, bands_plotted(b), d - 1)))... % drug_p_val_index(d) - 1)))...
                    
                    plot_test = reshape(plot_test, size(plot_test, 1), no_channels*2);
                    
                    plot_test(plot_test == 0) = nan;
                    
                    side_vec = [ones(1, 3) zeros(1, 3)];
                    
                    % side_vec(:, sum(~isnan(plot_test)) == 0) = [];
                    %
                    % plot_test(:, sum(~isnan(plot_test)) == 0) = [];
                    
                    subplot(no_bands_plotted, no_drugs - 2, (b - 1)*(no_drugs - 2) + d - 1)
                    
                    add_stars(gca, (1:timestep_lengths(t))', plot_test(1:timestep_lengths(t), :), side_vec, [])
                    
                    % add_stars(gca, (1:no_BP_hr_periods)', plot_test(:, :, 2), 0, [])
                    
                    % subplot(no_bands_plotted, no_drugs - 2, (b - 1)*(no_drugs - 2) + d - 1)
                    %
                    % y_lims = ylim;
                    %
                    % y_range = diff(y_lims);
                    %
                    % test_multiplier = [ones(size(plot_test(:, :, 2)))*diag(y_lims(1) - [.05 .1 .15]*y_range),...
                    %     ones(size(plot_test(:, :, 1)))*diag(y_lims(2) + [.15 .1 .05]*y_range)]; % [nan nan nan 0.05 .1 .15]*med_range);
                    %
                    % plot((1:no_BP_hr_periods)', [plot_test(:, :, 2) plot_test(:, :, 1)].*test_multiplier)
                    %
                    % ylim([y_lims(1) - .2*y_range, y_lims(2) + .2*y_range])
                    
                end
                
                linkaxes(ax(b, :))
                
            end
            
            % linkaxes([flipud(ax(:, 1)); ax(:, 2)])
            
            save_as_pdf(gcf,['ALL_BP',norms{n},'_multichannel_MK801_NVP', timesteps{t}]) %, 'orientation', 'portrait')
            
        end
        
    end
    
end

end
