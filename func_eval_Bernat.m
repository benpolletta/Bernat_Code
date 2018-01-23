function func_eval_Bernat(channels, drugs, subjects, fcn_handle, varargin)

present_dir = pwd;

for c = 1:length(channels)
    
    channel = channels{c};
    
    for d = 1:length(drugs)
        
        drug = drugs{d};
        
        for s = 1:length(subjects)
            
            record_dir = [subjects{s}, '_', drug];
            
            cd (record_dir)
            
            epoch_list = [record_dir, '_chan1_epochs.list'];
            
            epoch_names = text_read(epoch_list, '%s');
            
            epoch_no = length(epoch_names);
            
            first_data = load(epoch_names{1});
            
            first_output = feval(fcn_handle, first_data, varargin{:});
            
            output_size = size(first_output);
            
            % Find size of array needed to hold outputs from all windows, based on size
            % of output from a single window and on the number of windows in each
            % dimension.
            
            all_output_size = [output_size, epoch_no];
            
            output_dim = length(output_size);
            
            all_output_dim = output_dim + 1;
            
            % Initialize cell holding indices for each window, and size arguments for
            % array to hold outputs from all windows.
            
            all_output = nan(all_output_size); % Initialize array to hold outputs from all windows.
            
            output_indices = cell(output_dim, 1);
            
            output_indices(:) = {':'};
            
            all_output = nan(all_output_size_cell{:});
            
            all_output(output_indices{:}, 1) = first;
            
            for e = 1:epoch_no
                
                data = load(epoch_names{e});
                
                all_output(output_indices{:}, e) = feval(fcn_handle, data, varargin{:});
                
            end
            
            cd (present_dir)
            
            save([record_dir, '_chan1_whm.mat'], 'all_output')
            
        end
        
    end
    
end
