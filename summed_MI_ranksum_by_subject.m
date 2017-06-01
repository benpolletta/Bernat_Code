function summed_MI_ranksum_by_subject

load('channels.mat'), no_channels = length(channel_names);

load('drugs.mat'), no_drugs = length(drugs);

load('subjects.mat'), no_subjects = length(subjects);

norms = {'_', '_pct_'}; no_norms = length(norms);

timesteps = {'hr', '6min'}; no_timesteps = length(timesteps);

mat_labels = {'by_state', 'BP_stats'};

All_BP_ranksum = nan(10*8, no_channels, 2, 6, no_drugs - 1, no_norms, no_timesteps);

%% Plots by time (not broken down by state).

for n = 1:no_norms
    
    for t = 1:no_timesteps
        
        suffix = ['hrMI', norms{n}, timesteps{t}, '_', mat_labels{2}];
        
        clear All_BP_stats
        
        for ch = 1:no_channels
            
            ch_dir = ['ALL_', channel_names{ch}];
            
            for s = 1:no_subjects
                
                subject = subjects{s};
        
                load([ch_dir, '/', ch_dir, '_', subject, '_summed_', suffix, '.mat'])
                
                BP_stats_new = permute(BP_stats, [2, 1, 3, 4]);
                
                BPs_dims = size(BP_stats_new);
                
                BP_stats_new = reshape(BP_stats_new, BPs_dims(1), 1, 1, BPs_dims(2), BPs_dims(3), BPs_dims(4));
                
                All_BP_stats(:, s, ch, :, :, :) = BP_stats_new;
                
            end
            
            for d = 2:no_drugs
                
                for b = 1:size(All_BP_stats, 4)
                    
                    for pd = 1:size(All_BP_stats, 1)
                        
                        comparison_data = [All_BP_stats(pd, :, ch, b, 1, 1)' All_BP_stats(pd, :, ch, b, 1, d)'];
                        
                        if sum(any(~isnan(comparison_data))) == 2
                            
                            All_BP_ranksum(pd, ch, 1, b, d - 1, n, t) = ranksum(comparison_data(:, 1), comparison_data(:, 2), 'tail', 'left');
                            
                            All_BP_ranksum(pd, ch, 2, b, d - 1, n, t) = ranksum(comparison_data(:, 1), comparison_data(:, 2), 'tail', 'right');
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

save('summed_MI_ranksum_by_subject.mat', 'All_BP_ranksum')

end
