function A103_A104_saline_NVP_chan1_5_replace_spec

channel_labels = {'Frontal', 'CA1'};

channels = [1 5];

subject_labels={'A103','A104'};
subj_num=length(subject_labels);

drug_labels={'saline','NVP'};
drug_num=length(drug_labels);

for d = 1:drug_num
    
    drug = char(drug_labels{d});

    for s = 1:subj_num
        
        subject = char(subject_labels{s});
        
        channel_label = channel_labels{s};
        
        dir=['ALL_',channel_label];
        
        drugs = text_read([dir,'/',dir,'_drugs.txt'], '%s');
        subjects = text_read([dir,'/',dir,'_subjects.txt'], '%s');
        % hrs_fid=fopen([dir,'/',dir,'_hrs.txt']);
        % fourhrs_fid=fopen([dir,'/',dir,'_4hrs.txt']);
        % states_fid=fopen([dir,'/',dir,'_states.txt']);
        spec = load([dir,'/',dir,'_spec.txt']);
        spec_pct = load([dir,'/',dir,'_spec_pct.txt']);
        spec_zs = load([dir,'/',dir,'_spec_zs.txt']);
        
        save([dir,'/',dir,'_spec_OLD_',datestr(now,'mm-dd-yy_HH-MM-SS'),'.mat'], 'spec', 'spec_pct', 'spec_zs')
        
        channel = channels(s);
        
        record_dir = [subject,'_',drug];
        
        channel_name = [record_dir,'_chan',num2str(channel)];
        
        all_dir=['ALL_',subject,'_chan',num2str(channel)];
        
        fourhrs = text_read([all_dir,'/ALL_',channel_name,'_states_pds.txt'],'%*s%*s%s%*[^\n]');
        
        spec_mat = load([all_dir,'/ALL_',channel_name,'_spec.mat']);
        spec_data = spec_mat.spec_all;
       
        data_indices = strcmp(subjects, subject) & strcmp(drugs, drug);
        size(spec), sum(data_indices), size(spec_data)
        
        spec_data = permute(spec_data, [3 2 1]);
        sd_dims = size(spec_data);
        spec_data = reshape(spec_data(:, 1:(end - 1), :), 2, (sd_dims(2) - 1)/2, sd_dims(3));
        spec_data = permute(mean(spec_data), [3 2 1]);
        size(spec_data)
        
        spec(data_indices, :) = spec_data;
        
        [epochs, no_freqs] = size(spec_data);
        spec_format = make_format(no_freqs,'f');
        
        [~, est_inj_epoch] = max(spec_data(:, 1));
        spec_data(est_inj_epoch-10:est_inj_epoch+10, :) = nan;
        
        spec_mean = ones(epochs, no_freqs)*diag(nanmean(spec_data));
        spec_std = ones(epochs, no_freqs)*diag(nanmean(spec_data));
        spec_data_zs = (spec_data - spec_mean)./spec_std;
        spec_zs(data_indices, :) = spec_data_zs;
        
        baseline_indices = strcmp(fourhrs, 'pre4to1');
        baseline_spec = ones(epochs, no_freqs)*diag(nanmean(spec_data(baseline_indices, :)));
        spec_data_pct = 100*spec_data./baseline_spec - 100*ones(epochs, no_freqs);
        spec_pct(data_indices, :) = spec_data_pct;
        
    end
    
end

norms = {'', '_pct', '_zs'};

fid = nan(1, length(norms));

for n = 1:length(norms)

    fid(n) = fopen([dir,'/',dir,'_spec', norms{n}, '.txt'], 'w');

end
    
fprintf(fid(1), spec_format, spec');
fprintf(fid(2), spec_format, spec_pct');
fprintf(fid(3), spec_format, spec_zs');

fclose('all');
