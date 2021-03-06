function run_peak_avg_wt_high_MI_hrs_by_subject(start_index)

home_dir = '/projectnb/crc-nak/brpp';

[peak_hours, max_phase_freq, ~] = hour_peak_summed_MI_by_subject;

peak_freq_cycles = 3; sampling_freq = 1000;

load('subjects')
load('drugs')
load('channels')

if isempty(start_index), start_index = 0; end

index = 1;

for run = 1:2
    
    for s = 1:no_subjects
        
        for d = 1:no_drugs
            
            hr_label = sprintf('post%d', peak_hours(s, d, run + 4));
            
            hr_phase = max_phase_freq(s, d, run + 4);
            
            record_dir = [subjects{s}, '_', drugs{d}, '_', hr_label, '_epochs'];
            
            cd ([home_dir, '/Bernat_NMDAR_antagonists/', record_dir])
            
            load([record_dir, '_numbers.mat'])
            
            for c = 1:no_channels
                
                if index > start_index
                    
                    channel = location_channels{c}(s);
                    
                    channel_listname = [record_dir, '_chan', num2str(channel), '.list'];
                    
                    fid = fopen(channel_listname, 'w');
                    
                    for e = 1:length(epoch_numbers)
                        
                        epoch_filename = [subjects{s}, '_', drugs{d}, '_chan', num2str(channel), '_epoch', num2str(epoch_numbers(e)), '.txt'];
                        
                        fprintf(fid, '%s\n', epoch_filename);
                        
                    end
                    
                    fclose(fid);
                    
                    peak_averaged_wt_batch_parallel(channel_listname, hr_phase, peak_freq_cycles, sampling_freq);
                    
                end
                
                index = index + 1;
                
            end
                    
            cd ([home_dir, '/Bernat_NMDAR_antagonists'])
            
        end
        
    end
    
end

end

% function max_phase = get_peak_freqs
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
% All_cplot_data = nan(no_afs, no_pfs, no_drugs, no_hr_periods, no_channels);
% 
% for c = 1:no_channels
%     
%     for d = 1:no_drugs
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
% max_phase_data = reshape(max(All_cplot_data), no_pfs, no_drugs, no_hr_periods, no_channels);
% 
% max_phase = nan(no_hr_periods, no_drugs, no_channels);
% 
% for c = 1:no_channels
%     
%     for d = 1:no_drugs
%         
%         for h = 1:no_hr_periods
%            
%             phase_profile = max_phase_data(:, d, h, c);
%             
%             [~, max_index] = max(phase_profile);
%             
%             max_phase(h, d, c) = phase_freqs(max_index);
%             
%         end
%         
%     end
%     
% end
% 
% % max_phase = reshape(median(max_phase), no_drugs, no_hr_periods);
% 
% end