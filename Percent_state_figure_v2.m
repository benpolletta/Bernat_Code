function Percent_state_figure_v2(significance)

if nargin < 1, significance = []; end
if isempty(significance), significance = .025; end
sig_flag = make_label('p', significance, .025);

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat'), no_drugs = length(drugs);

load('states.mat'), no_states = length(states);
long_state_labels = {'Active Wake', 'Quiet Wake/nREM', 'REM'};

stats={'mean'}; % dian'}; %,'mean','std'}; 
no_stats = length(stats);
long_stats={'Mean'}; % ,'Mean','St. Dev.'};

timesteps = {'_hrs', '_6mins'};
no_timesteps = length(timesteps);

no_pre=4; no_post=20;
[BP_hr_labels, ~, ~] = make_period_labels(no_pre, no_post, 'hrs');
no_BP_hr_periods = length(BP_hr_labels);
short_BP_hr_labels = -no_pre:no_post; short_BP_hr_labels(short_BP_hr_labels == 0) = [];

no_pre=2; no_post=8;
[BP_6min_labels, ~, ~] = make_period_labels(no_pre, no_post, '6mins');
no_BP_6min_periods = length(BP_6min_labels);
short_BP_6min_labels = ((-no_pre*60 + 3):6:(no_post*60 - 3))/60; short_BP_6min_labels(short_BP_6min_labels == 0) = [];

timestep_labels = {short_BP_hr_labels, short_BP_6min_labels};
timestep_lengths = cellfun(@(x) length(x), timestep_labels);
tick_spacing = [4 2]; % floor(timestep_lengths/4);

drug_p_val_index = [1 4 2 3];
    
clear titles xlabels ylabels

load('Percent_state_ranksum_by_subject')

last_drug = no_drugs; % - 1;

no_drugs_plotted = last_drug - 2 + 1;

% c_order = [.25 .25 .25; .5 .25 .25; 0 .5 0; .75 0 .75];
  
c_order = distinguishable_colors(11); c_order = [.25 .25 .25; c_order(9:11, :)];

for t = 1:no_timesteps
    
    clear All_BP_stats
    
    load(['Percent_state', timesteps{t}, '_BP_stats.mat'])
    
    % Testing p-values.
    
    All_BP_test = All_BP_ranksum(:, :, :, :, t) <= significance; % /(3*sum_all_dimensions(~isnan(All_BP_ranksum(:, :, :, :, :, n, t))));
    
    figure
    
    %% Plotting time series.
    
    for state = 1:no_states
        
        clear plot_mean plot_std
        
        % plot_median = BP_stats(:, :, 1, d) - BP_stats(:, :, 1, 1);
        
        plot_mean = squeeze(BP_stats(state, :, 4, :));
        
        plot_std = squeeze(BP_stats(state, :, 5, :));
        
        subplot(no_states, 1, state)
        
        set(gca, 'NextPlot', 'add', 'ColorOrder', c_order)
        
        plot(timestep_labels{t}', plot_mean*100, 'LineWidth', 2)
        
        % boundedline(timestep_labels{t}', plot_mean, prep_for_boundedline(plot_std), 'cmap', c_order)
        
        hold on
        
        set(gca, 'XTick', min(timestep_labels{t}):tick_spacing(t):max(timestep_labels{t}),...
            'XTickLabel', round(min(timestep_labels{t}):tick_spacing(t):max(timestep_labels{t})),...
            'FontSize', 16)
        
        axis tight
        
        ylabel(long_state_labels{state})
        
        if state == no_states
            
            xlabel('Time Rel. Inj. (h)')
            
        end
        
        if state == 1
            
            title('Percentage of Time')
            
            legend(drugs)
        
        end
        
        %% Plotting stats.
        
        clear plot_test
        
        plot_test = squeeze(All_BP_test(:, :, state, :)); %]; % drug_p_val_index(d) - 1)]; [nan(size(All_BP_test(:, :, bands_plotted(b), d - 1)))... % drug_p_val_index(d) - 1)))...
        
        plot_test = double(permute(plot_test, [1 3 2]));
        
        no_test_vars = size(plot_test, 2);
        
        % set(gca, 'NextPlot', 'add', 'ColorOrder', c_order(2:end, :))
        
        plot_test = reshape(plot_test, size(plot_test, 1), size(plot_test, 3)*no_test_vars);
        
        plot_test(plot_test == 0) = nan;
        
        side_vec = [ones(1, no_test_vars) zeros(1, no_test_vars)];
        
        % subplot(no_states, 1), state)
        
        add_stars(gca, timestep_labels{t}', plot_test(1:length(timestep_labels{t}), :), side_vec, c_order(2:end, :))
        
    end
    
    save_as_pdf(gcf,['Percent_state', timesteps{t}, '_figure_v2', sig_flag]) %, 'orientation', 'portrait')
    
end
