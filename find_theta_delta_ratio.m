function find_theta_delta_ratio(drug)

load('subjects.mat'), load('channels.mat')

ch_selected = [1 3];

freq_ranges = {[0.25 4.75], [5.25 10.75]};

present_dir = pwd;

for s = 1:length(subjects)
    
    subject = subjects{s};
    
    record_dir = [subjects{s}, '_', drug];
    
    cd (record_dir)
    
    for ch = 1:2
        
        ch_i = ch_selected(ch);
        
        ch_no = location_channels{ch_i}(s);
        
        epoch_list = sprintf('%s_chan%d_epochs.list', record_dir, ch_no);
        
        epoch_names{ch} = text_read(epoch_list, '%s');
        
        no_epochs(ch) = length(epoch_names{ch});
        
    end
    
    if diff(no_epochs) ~= 0
        
        fprintf('For %s, the number of epochs for channel %d and channel %d are different (%d vs. %d).\n',...
            subjects{s}, location_channels{ch_selected(1)}(s), location_channels{ch_selected(2)}(s), no_epochs)
        
        cd (present_dir)
       
        return
        
    end
    
    [power, max_freqs, max_vals] = deal(nan(no_epochs(1), 2));
    
    parfor e = 1:no_epochs(1)
        
        for ch = 1:2
        
            data = load(epoch_names{ch}{e});
            
            [~, ~, max_freqs(e, ch), max_vals(e, ch), ~, power(e, ch)] = width_half_max(data, 1000, freq_ranges{ch}, .075, 0)
            
        end
        
    end
    
    val_ratio = max_vals(:, 1)./max_vals(:, 2);
    
    pow_ratio = power(:, 1)./power(:, 2);
    
    cd (present_dir)
    
    save([record_dir, '_tdr.mat'], 'max_freqs', 'max_vals', 'power')
    
    tdr_mat = [val_ratio pow_ratio power]; no_criteria = size(tdr_mat, 2);
    
    crit_labels = {'val_ratio', 'pow_ratio', 'delta', 'theta'};
    
    figure
    
    for c = 1:no_criteria
        
        subplot(2, no_criteria, c)
        
        [n, centers] = hist(tdr_mat(:, c), 25);
        
        plot(centers, n)
        
        axis tight
        
        title(crit_labels{c})
        
        subplot(2, no_criteria, no_criteria + c)
        
        plot(tdr_mat(:, c), '*')
        
        axis tight
        
    end
    
    save_as_pdf(gcf, [record_dir, '_tdr'])
    
end