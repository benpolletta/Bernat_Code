function plot_entropy_time_series(drug, quantile_used, states)

load('subjects.mat'), load('AP_freqs.mat')
    
state_label = '';

if ~isempty(states)
    
    no_states = length(states);
    
    for state = 1:no_states
        
        state_label = [state_label, '_', states{state}];
        
    end
    
end

name = 'ALL_Frontal';

measure = 'p0.99_IEzs';

MI_drugs = text_read([name,'/',name,'_',measure,'_drugs.txt'],'%s');
MI_subjects = text_read([name,'/',name,'_',measure,'_subjects.txt'],'%s');
MI_6mins = text_read([name, '/', name, '_', measure, '_6mins.txt'], '%d');
MI_states = text_read([name,'/',name,'_',measure,'_states.txt'],'%s');

no_pre=4; no_post=8;
short_BP_6min_labels = ((-no_pre*60 + 3):6:(no_post*60 - 3));

for s = 1:subj_num
    
    subject = subjects{s};
    
    record_dir = [subject, '_', drug];
    
    subj_MI_index = strcmp(MI_subjects, subject) & strcmp(MI_drugs, drug);
    
    if ~isempty(states)
       
        subj_state_index = zeros(sum(subj_MI_index), 1);
        
        for state = 1:no_states
            
            subj_state_index = subj_state_index | strcmp(MI_states(subj_MI_index), states{state});
            
        end
        
    else
        
        subj_state_index = ones(sum(subj_MI_index), 1);
        
    end
    
    clear entropy
    
    load([record_dir, '_chan1_whm.mat'])
    
    if length(whm) > length(subj_state_index)
        
        entropy((length(subj_state_index) + 1):end) = [];
        
    end
    
    low_entropy_indices = entropy < quantile(entropy, quantile_used) & subj_state_index;
    
    high_entropy_indices = entropy > quantile(entropy, 1 - quantile_used) & subj_state_index;
    
    for sixmin = 1:length(short_BP_6min_labels)
        
        sixmin_indices = MI_6mins(subj_MI_index) == short_BP_6min_labels(sixmin);
        
        entropy_hist(sixmin, s, 1) = sum(low_entropy_indices(sixmin_indices));
        
        entropy_hist(sixmin, s, 2) = sum(high_entropy_indices(sixmin_indices));
        
    end
    
end

figure

subplot(3, 1, 1)

plot(short_BP_6min_labels/60, entropy_hist(:, :, 1))

title(drug, 'FontSize', 16)

ylabel('Low Entropy')

axis tight

subplot(3, 1, 2)

plot(short_BP_6min_labels/60, entropy_hist(:, :, 2))

ylabel('High Entropy')

axis tight

subplot(3, 1, 3)

plot(short_BP_6min_labels/60, squeeze(sum(entropy_hist, 2)))

axis tight

xlabel('Time Rel. Infusion (h)')

legend({'Low Entropy', 'High Entropy'})

save_as_pdf(gcf, sprintf('%s_entropy_q%g_%s', drug, quantile_used, state_label))

