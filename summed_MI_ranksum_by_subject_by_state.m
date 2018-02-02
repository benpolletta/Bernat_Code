function summed_MI_ranksum_by_subject_by_state

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat'), no_drugs = length(drugs);

load('subjects.mat'), no_subjects = length(subjects);

load('states.mat'), no_states = length(states);

norms = {'_', '_pct_'}; no_norms = length(norms);

no_pre=4; no_post=20;
[BP_hr_labels, ~, ~] = make_period_labels(no_pre, no_post, 'hrs');
no_BP_hr_periods = length(BP_hr_labels);
short_BP_hr_labels = -no_pre:no_post; short_BP_hr_labels(short_BP_hr_labels == 0) = [];

no_pre=2; no_post=8;
[BP_6min_labels, ~, ~] = make_period_labels(no_pre, no_post, '6mins');
no_BP_6min_periods = length(BP_6min_labels);
short_BP_6min_labels = ((-no_pre*60 + 3):6:(no_post*60 - 3))/60; short_BP_6min_labels(short_BP_6min_labels == 0) = [];

timesteps = {'hr', '6min'};
timestep_labels = {short_BP_hr_labels, short_BP_6min_labels};
no_timesteps = length(timestep_labels);
timestep_lengths = cellfun(@(x) length(x), timestep_labels);
% tick_spacing = floor(timestep_lengths/4);

mat_labels = {'by_state', 'BP_stats'};

no_bands = 6;

All_BP_stats_by_subject = nan(10*10, no_subjects, no_bands, no_states, 5, no_drugs,...
    no_norms, no_timesteps, no_channels);

All_BP_ranksum = nan(10*10, no_channels, 2, no_bands, no_drugs - 1,...
    no_norms, no_timesteps, no_states);

%% Plots by time (not broken down by state).

for n = 1:no_norms
    
    for t = 1:no_timesteps
        
        suffix = ['hrMI', norms{n}, timesteps{t}, '_', mat_labels{1}];
        
        clear All_BP_stats
        
        for ch = 1:no_channels
            
            ch_dir = ['ALL_', channel_names{ch}];
            
            for s = 1:no_subjects
                
                subject = subjects{s};
        
                load([ch_dir, '/', ch_dir, '_', subject, '_summed_', suffix, '.mat'])
                
                BP_stats_new = permute(BP_stats, [3, 1, 2, 4, 5]);
                
                BPs_dims = size(BP_stats_new);
                
                BP_stats_new = reshape(BP_stats_new, [BPs_dims(1), 1, BPs_dims(2:end)]);
                
                All_BP_stats(:, s, :, :, :, :) = BP_stats_new;
                
            end
            
            All_BP_stats_by_subject(1:timestep_lengths(t), :, :, :, :, :, n, t, ch) = All_BP_stats;
            
            for d = 2:no_drugs
                
                for b = 1:no_bands
                    
                    for state = 1:no_states
                        
                        for pd = 1:size(All_BP_stats, 1)
                            
                            comparison_data = [All_BP_stats(pd, :, b, state, 1, 1)' All_BP_stats(pd, :, b, state, 1, d)'];
                            
                            if sum(any(~isnan(comparison_data))) == 2
                                
                                All_BP_ranksum(pd, ch, 1, b, d - 1, n, t, state) = ranksum(comparison_data(:, 1), comparison_data(:, 2), 'tail', 'left');
                                
                                All_BP_ranksum(pd, ch, 2, b, d - 1, n, t, state) = ranksum(comparison_data(:, 1), comparison_data(:, 2), 'tail', 'right');
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

save('summed_MI_ranksum_by_subject_by_state.mat', 'All_BP_ranksum', 'All_BP_stats_by_subject')

end
