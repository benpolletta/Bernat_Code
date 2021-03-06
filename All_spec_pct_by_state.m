function All_spec_pct_by_state(channel_label)

subject_labels = load('subjects');
subject_labels = subject_labels.subjects;
subj_num = length(subject_labels);

drug_labels = load('drugs');
drug_labels = drug_labels.drugs;
drug_num = length(drug_labels);

state_labels = load('states');
state_labels = state_labels.states;
no_states = length(state_labels);

dir = ['ALL_',channel_label];

drugs = text_read([dir,'/',dir,'_drugs.txt'], '%s');
subjects = text_read([dir,'/',dir,'_subjects.txt'], '%s');
fourhrs = text_read([dir,'/',dir,'_4hrs.txt'], '%s');
states = text_read([dir,'/',dir,'_states.txt'], '%s');
spec = load([dir,'/',dir,'_spec.txt'],'w');
BP = load([dir,'/',dir,'_BP.txt'],'w');

spec_pct = nan(size(spec));
BP_pct = nan(size(BP));

for st = 1:no_states
    
    state = char(state_labels{st});
    
    state_indices = strcmp(states, state);
    
    for d = 1:drug_num
        
        drug = char(drug_labels{d});
        
        drug_indices = strcmp(drugs, drug);
        
        for s = 1:subj_num
            
            subject = char(subject_labels{s});
            
            subj_indices = strcmp(subjects, subject) & drug_indices & state_indices;
            
            subj_spec_pct = spec(subj_indices, :);
            
            subj_baseline_indices = strcmp(fourhrs, 'pre4to1') & subj_indices;
            
            baseline_spec = ones(size(subj_spec_pct))*diag(nanmean(spec(subj_baseline_indices,:)));
            
            subj_spec_pct = 100*subj_spec_pct./baseline_spec - 100*ones(size(subj_spec_pct));
            
            spec_pct(subj_indices, :) = subj_spec_pct;
            
            subj_BP_pct = BP(subj_indices, :);
            
            baseline_BP = ones(size(subj_BP_pct))*diag(nanmean(BP(subj_baseline_indices,:)));
            
            subj_BP_pct = 100*subj_BP_pct./baseline_BP - 100*ones(size(subj_BP_pct));
            
            BP_pct(subj_indices, :) = subj_BP_pct;
            
        end
        
    end
    
end

save([dir,'/',dir,'_spec_pct_by_state.mat'], 'spec_pct', 'BP_pct')