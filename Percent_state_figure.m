function Percent_state_figure(significance)

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
tick_spacing = floor(timestep_lengths/4);

drug_p_val_index = [1 4 2 3];

c_order = [0 0 1; 0 .5 0; 1 0 0];
    
clear titles xlabels ylabels

load('Percent_state_ranksum_by_subject')

last_drug = no_drugs; % - 1;

no_drugs_plotted = last_drug - 2 + 1;
    
for t = 1:no_timesteps
    
    clear All_BP_stats
    
    load(['Percent_state', timesteps{t}, '_BP_stats.mat'])
    
    % Testing p-values.
    
    All_BP_test = All_BP_ranksum(:, :, :, :, t) <= significance; % /(3*sum_all_dimensions(~isnan(All_BP_ranksum(:, :, :, :, :, n, t))));
    
    figure
    
    %% Plotting time series.
    
    for state = 1:no_states
        
        for d = 2:last_drug
            
            clear plot_stats
            
            % plot_median = BP_stats(:, :, 1, d) - BP_stats(:, :, 1, 1);

            plot_mean = squeeze(BP_stats(state, :, 4, [1 d]));
            
            plot_std = squeeze(BP_stats(state, :, 5, [1 d]));
            
            ax(state, d - 1) = subplot(no_states, no_drugs_plotted, (state - 1)*no_drugs_plotted + d - 1);
            
            % set(gca, 'NextPlot', 'add', 'LineStyleOrder', {'-','*','*'}, 'ColorOrder', c_order)
            %
            % plot((1:timestep_lengths(t))', plot_stats(1:timestep_lengths(t), :))
            
            boundedline(timestep_labels{t}', plot_mean, prep_for_boundedline(plot_std))
            
            hold on
            
            set(gca, 'XTick', min(timestep_labels{t}):tick_spacing(t):max(timestep_labels{t}),...
                'FontSize', 16)
            
            axis tight
            
            if state == 1
                
                title(drugs{d})
                
            elseif state == no_states
                
                xlabel('Time Rel. Inj. (h)')
                
            end
            
            if d - 1 == 1
                
                ylabel(long_state_labels{state})
                
            end
            
        end
        
        sync_axes(ax(state, :))
        
        %% Plotting stats.
        
        for d = 2:last_drug
            
            clear plot_test
            
            plot_test = double(All_BP_test(:, :, state, d - 1)); %]; % drug_p_val_index(d) - 1)]; [nan(size(All_BP_test(:, :, bands_plotted(b), d - 1)))... % drug_p_val_index(d) - 1)))...
            
            % plot_test = double(reshape(plot_test, size(plot_test, 1), no_channels*2));
            
            plot_test(plot_test == 0) = nan;
            
            side_vec = [1 0]; % [ones(1, no_channels) zeros(1, no_channels)];
            
            subplot(no_states, no_drugs_plotted, (state - 1)*no_drugs_plotted + d - 1)
            
            add_stars(gca, timestep_labels{t}', plot_test(1:length(timestep_labels{t}), :), side_vec, [])
            
        end
        
        sync_axes(ax(state, :))
        
    end
    
    save_as_pdf(gcf,['Percent_state', timesteps{t}, '_figure', sig_flag]) %, 'orientation', 'portrait')
    
end
