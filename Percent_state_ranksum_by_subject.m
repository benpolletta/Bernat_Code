function Percent_state_ranksum_by_subject

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat'), no_drugs = length(drugs);

load('subjects.mat'), no_subjects = length(subjects);

load('states.mat'), no_states = length(states);

timesteps = {'hrs', '6mins'}; no_timesteps = length(timesteps);

All_BP_ranksum = nan(10*8, 2, no_states, no_drugs - 1, no_timesteps);

for t = 1:no_timesteps
    
    suffix = [timesteps{t}, '_BP_stats'];
    
    clear All_BP_stats
    
    %% Loading percent state (mean over state indicator for 6 minute period) for each subject.
    
    for subj = 1:no_subjects
        
        subject = subjects{subj};
        
        load(['Percent_state_', subject, '_', suffix, '.mat']) % Dimensions 3x100x5x4
        
        BP_stats_new = permute(BP_stats, [2, 1, 3, 4]); % Now 100x3x5x4
        
        BPs_dims = size(BP_stats_new);
        
        BP_stats_new = reshape(BP_stats_new, BPs_dims(1), 1, BPs_dims(2), BPs_dims(3), BPs_dims(4)); % Now 100x1x1x3x5x4
        
        All_BP_stats(:, subj, :, :, :) = BP_stats_new; % Now 100x6x3x5x4
        
    end
    
    %% Comparing drugs to saline, by period, using percent time in each state as an observation for each subject.
    
    for d = 2:no_drugs
        
        for state = 1:size(All_BP_stats, 3)
            
            for pd = 1:size(All_BP_stats, 1)
                
                comparison_data = [All_BP_stats(pd, :, state, 4, 1)' All_BP_stats(pd, :, state, 4, d)'];
                
                if sum(any(~isnan(comparison_data))) == 2
                    
                    All_BP_ranksum(pd, 1, state, d - 1, t) = ranksum(comparison_data(:, 1), comparison_data(:, 2), 'tail', 'left');
                    
                    All_BP_ranksum(pd, 2, state, d - 1, t) = ranksum(comparison_data(:, 1), comparison_data(:, 2), 'tail', 'right');
                    
                end
                
            end
            
        end
        
    end
    
end

save('Percent_state_ranksum_by_subject.mat', 'All_BP_ranksum')

end
