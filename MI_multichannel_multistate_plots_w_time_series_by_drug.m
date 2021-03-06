function MI_multichannel_multistate_plots_w_time_series_by_drug(hi_hr, cplot_norm, drug)

phase_freqs = 1:.25:12; amp_freqs = 20:5:200;
no_pfs = length(phase_freqs); no_afs = length(amp_freqs);

load('channels.mat'), no_channels = length(channel_names);

load('states.mat'), no_states = length(states);

load('drugs.mat')

stats={'median','mean','std'}; no_stats = length(stats);
long_stats={'Median','Mean','St. Dev.'}; total_stats = 5;

norms={'', '_pct'}; no_norms = length(norms);
long_norms={'', '% Change'};% From Baseline'};

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

c_order = [0 0 1; 0 .5 0; 1 0 0];
    
clear titles xlabels ylabels

preinj_data = nan(no_afs, no_pfs, no_stats, no_channels, no_norms);
           
All_cplot_data = nan(no_afs, no_pfs, no_states, no_hr_periods, no_stats, no_channels, no_norms);

All_BP_stats = nan(no_BP_hr_periods, no_states, no_channels, no_bands, total_stats, no_norms);

[All_BP_ranksum, All_BP_test] = deal(nan(no_BP_hr_periods, no_states - 1, no_channels, no_bands, no_norms));
    
for n=1:no_norms
    
    for c=1:no_channels
        
        ch_dir = ['ALL_', channel_names{c}];
            
        %% Getting colorplot data.
            
        load([ch_dir,'/',ch_dir,'_p0.99_IEzs_MI','/',...
            ch_dir,'_p0.99_IEzs_hrMI_hr', norms{n},'_by_state_',drug,'_cplot_data.mat'])
        
        All_cplot_data(:, :, :, :, :, c, n) = permute(MI_stats, [1 2 3 4 6 5]);
        
        %% Getting time series data.
        
        suffix = ['4hrMI', norms{n}, '_hr_by_state'];
        
        load([ch_dir, '/', ch_dir, '_summed_', suffix, '.mat'])
        
        BP_drug_stats = BP_stats(:, :, 1:no_BP_hr_periods, :, strcmp(cat3_labels, drug));
        
        All_BP_stats(:, :, c, :, :, n) = permute_fit(All_BP_stats(:, :, c, :, :, n), BP_drug_stats);
        
        %% Getting stats p-values.
        
        if strcmp(norms{n}, '_pct')
        
            ranksum_suffix = 'hrMI_pct_by_state_hr_state_vs_wake_ranksum';
        
        else
            
            ranksum_suffix = 'hrMI_hr_state_vs_wake_ranksum';
            
        end
        
        load([ch_dir, '/', ch_dir, '_summed_', ranksum_suffix, '.mat'])
        
        BP_drug_ranksum = BP_ranksum(:, strcmp(cat3_labels, drug), :, 1:no_BP_hr_periods);
        
        All_BP_ranksum(:, :, c, :, n) = permute_fit(All_BP_ranksum(:, :, c, :, n), BP_drug_ranksum);
        
    end
    
    % Bonferroni correcting & testing p-values.
    
    All_BP_test(:, :, :, :, n) = All_BP_ranksum(:, :, :, :, n) <= .01/(3*sum_all_dimensions(~isnan(All_BP_ranksum(:, :, :, :, n))));
    
end
    
All_BP_stats(:, :, :, :, [2 3], :) = []; % Leaving behind [1 4 5], which are median, mean, std.

if strcmp(hi_hr, 'independent')

    [~, max_hr_indices] = nanmax(nanmax(nanmax(All_cplot_data(:, :, :, 4:end, :, :, :))), [], 4);
    
    max_hr_indices = reshape(max_hr_indices, no_states, no_stats, no_channels, no_norms) + 4 - 1;

elseif strcmp(hi_hr, 'state')
   
    [~, max_hr_indices] = nanmax(nanmax(nanmax(nanmax(All_cplot_data(:, :, :, 4:end, :, :, :))), [], 6), [], 4);
    
    max_hr_indices = repmat(reshape(max_hr_indices, no_states, no_stats, 1, no_norms), [1 1 no_channels 1]) + 4 - 1;
    
end

handle = nan(no_states, no_norms, no_stats);

All_cplot_for_plot = nan(no_afs, no_pfs, no_states, no_stats, no_channels, no_norms);

for n = 1:no_norms
    
    for c = 1:no_channels
        
        for s = 1:no_stats
            
            for st = 1:no_states
                
                All_cplot_for_plot(:, :, st, s, c, n) = All_cplot_data(:, :, st, max_hr_indices(st, s, c, n), s, c, n);
                
            end
            
        end
        
    end
    
end

if strcmp(cplot_norm, '_row')
    
    max_by_channel = reshape(nanmax(nanmax(nanmax(All_cplot_for_plot))), no_stats, no_channels, no_norms);
    
    min_by_channel = reshape(nanmin(nanmin(nanmin(All_cplot_for_plot))), no_stats, no_channels, no_norms);
    
elseif strcmp(cplot_norm, '_col')
    
    max_by_drug = reshape(nanmax(nanmax(nanmax(All_cplot_for_plot)), [], 5), no_states, no_stats, no_norms);
    
    min_by_drug = reshape(nanmin(nanmin(nanmin(All_cplot_for_plot)), [], 5), no_states, no_stats, no_norms);
    
end

for n = 1:no_norms
    
    for s = 1:no_stats
        
        handle(n, s) = figure;
        
        for c = 1:no_channels
            
            for st = 1:no_states
                
                %% Plotting comodulograms.
                
                subplot(no_channels + no_bands_plotted, no_states, (c - 1)*no_states + st)
                
                imagesc(phase_freqs, amp_freqs, All_cplot_for_plot(:, :, st, s, c, n))
                
                if strcmp(cplot_norm, '_row')
                    
                    caxis([min_by_channel(s, c, n) max_by_channel(s, c, n)])
                    
                    if st == no_states
                        
                        colorbar
                        
                    end
                    
                elseif strcmp(cplot_norm, '_col')
                    
                    caxis([min_by_drug(st, s, n) max_by_drug(st, s, n)])
                    
                else
                    
                    colorbar
                    
                end
                
                axis xy
                
                if c == 1
                    
                    title({[states{st}, ', ', drug]; [long_stats{s}, ' MI, ', long_norms{n}]; long_hr_labels{max_hr_indices(st, s, c, n)}})
                    
                else
                    
                    title(long_hr_labels{max_hr_indices(st, s, c, n)})
                    
                end
                
                if st == 1
                    
                    ylabel(channel_names{c})
                    
                end
                
            end
            
            %% Plotting time series w/ stats.
            
            for b = 1:no_bands_plotted
                
                subplot(2*no_channels, no_bands_plotted, no_channels*no_bands_plotted + (c - 1)*no_bands_plotted + b)
                
                clear plot_stats plot_test
                
                plot_stats = All_BP_stats(:, :, c, bands_plotted(b), s, n);
                
                % plot((1:no_BP_hr_periods)', plot_stats)
                
                % add_stars(gca, (1:no_BP_hr_periods)', logical(All_BP_test(:, :, c, bands_plotted(b), n)), [1 1], c_order);
                
                plot_test = [nan(size(All_BP_test, 1), 1) All_BP_test(:, :, c, bands_plotted(b), n)];
                
                plot_test(plot_test == 0) = nan;
                
                med_min = min(min(plot_stats(:, :, 1)));
                
                med_range = max(max(plot_stats(:, :, 1))) - med_min;
                
                test_multiplier = ones(size(plot_test))*diag(med_min - [nan 0.05 .1]*med_range);
                
                % subplot(no_channels + no_bands_plotted, no_states - 1, (no_channels + (b - 1))*(no_states - 1) + d - 1)
                
                set(gca, 'NextPlot', 'add', 'LineStyleOrder', {'-','*'}, 'ColorOrder', c_order)
                
                plot((1:no_BP_hr_periods)', [plot_stats plot_test.*test_multiplier])
                
                set(gca, 'XTick', 1:tick_spacing:no_BP_hr_periods, 'XTickLabel', short_BP_hr_labels(1:tick_spacing:end))
                
                axis tight
                
                ylim([med_min - .2*med_range, med_min + 1.05*med_range])
                
                if c == 1
                    
                    title(band_labels{bands_plotted(b)})
                    
                    if b == 1
                        
                        legend(states, 'Location', 'NorthEast', 'FontSize', 6)
                        
                    end
                    
                elseif c == no_channels
                    
                    xlabel('Time Rel. Inj. (h)')
                    
                end
                
                if b == 1
                    
                    ylabel(channel_names{c})
                    
                end
                
            end
            
        end
                
        save_as_pdf(gcf,['ALL_MI',norms{n},'_multichannel_multistate_', hi_hr, '_hi_', drug, '_', stats{s}, cplot_norm])
        
    end
    
end
    
% for n=1:no_norms
%     
%     for s=1:no_stats
%             
%         open(['ALL_MI',norms{n},'_multichannel_multistate_', hi_hr, '_hi_', drug, '_', stats{s}, cplot_norm, '.fig'])
%         
%     end
%     
% end

end

% function [max_phase, max_pd] = get_peak_periods
% 
% phase_freqs = 1:.25:12; amp_freqs = 20:5:200;
% no_pfs = length(phase_freqs); no_afs = length(amp_freqs);
% 
% load('channels.mat'), no_channels = length(channel_names);
% 
% load('drugs.mat')
% 
% stats={'median','mean','std'}; no_stats = length(stats);
% 
% no_pre=4; no_post=12;
% [hr_labels, ~, ~]=make_period_labels(no_pre,no_post,'hrs');
% no_hr_periods = length(hr_labels);
% 
% All_cplot_data = nan(no_afs, no_pfs, no_states, no_hr_periods, no_channels);
% 
% for c = 1:no_channels
%     
%     for st = 1:no_states
%         
%         load(['ALL_',channel_names{c},'/ALL_',channel_names{c},'_p0.99_IEzs_MI',...
%             '/ALL_',channel_names{c},'_p0.99_IEzs_hr_',drugs{d},'_cplot_data.mat'])
%         
%         All_cplot_data(:, :, d, :, c) = MI_stats(:, :, 1, :, 1);
%         
%     end
%     
% end
% 
% max_MI_data = reshape(max(max(All_cplot_data)), no_states, no_hr_periods, no_channels);
% 
% [~, max_pd_indices] = max(max_MI_data, [], 2);
% 
% max_pd = hr_labels(reshape(max_pd_indices, no_states, no_channels));
% 
% max_phase_data = reshape(max(All_cplot_data), no_pfs, no_states, no_hr_periods, no_channels);
%  
% [~, max_phase_indices] = max(max_phase_data);
% 
% max_phase = permute(phase_freqs(max_phase_indices), [2 1 3]);
% 
% end
