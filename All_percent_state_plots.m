function All_percent_state_plots

subject_labels={'A99','A102','A103','A104','A105','A106'};
subj_num=length(subject_labels);

state_labels = {'W','NR','R'};
long_state_labels = {'Active Wake', 'Quiet Wake/nREM', 'REM'};
no_states = length(state_labels);

drug_labels = {'saline','MK801','NVP','Ro25'};
no_drugs = length(drug_labels);

name = 'ALL_Frontal';

drugs = text_read(sprintf('%s/%s_drugs.txt', name, name), '%s');
subjects = text_read(sprintf('%s/%s_subjects.txt', name, name), '%s');
states = text_read(sprintf('%s/%s_states.txt', name, name), '%s');

state_indicator = nan(length(states), length(state_labels));

for s = 1:no_states

    state_indicator(:, s) = cellfun(@(x) strcmp(x, state_labels{s}), states);
    
end

c_order = [.5 .5 .5; 1 .75 0; 1 0 1; 0 1 1];

timesteps = {'hrs', '6mins'};
no_timesteps = length(timesteps);

no_pre = [4 2]; no_post = [20 8];

for t = 1:no_timesteps
    
    pds = text_read(sprintf('%s/%s_%s.txt', name, name, timesteps{t}), '%s');
    
    [pd_labels, ~] = make_period_labels(no_pre(t), no_post(t), timesteps{t});
    
    lineplot_collected_BP_by_categories('Vigilance State',...
        sprintf('Percent_state_%s', timesteps{t}),...
        {state_labels, long_state_labels},{drug_labels, drug_labels}, {pd_labels, pd_labels},...
        drugs, pds, state_indicator, 'Mean', c_order)
    
    for subj = 1:subj_num
        
        subj_indicator = strcmp(subjects, subject_labels{subj});
        
        lineplot_collected_BP_by_categories([subject_labels{subj}, ', Vigilance State'],...
            sprintf('Percent_state_%s_%s', subject_labels{subj}, timesteps{t}),...
            {state_labels, long_state_labels},{drug_labels, drug_labels}, {pd_labels, pd_labels},...
            drugs(subj_indicator), pds(subj_indicator), state_indicator(subj_indicator, :), 'Mean', c_order)
        
    end
    
end