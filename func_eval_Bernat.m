function func_eval_Bernat(channels, drugs, subjects, fcn_handle, no_args_out, varargin)

present_dir = pwd;

load('channels.mat')

for c = 1:length(channels)
    
    channel_number = find(strcmp(channels{c}, channel_names));
    
    for d = 1:length(drugs)
        
        drug = drugs{d};
        
        for s = 1:length(subjects)
            
            record_dir = [subjects{s}, '_', drug];
            
            cd (record_dir)
            
            epoch_list = [record_dir, sprintf('_chan%d_epochs.list', location_channels{channel_number}(s))];
            
            epoch_names = text_read(epoch_list, '%s');
            
            epoch_no = length(epoch_names);
            
            first_data = load(epoch_names{1});
            
            [first_outputs{1:no_args_out}] = feval(fcn_handle, first_data, varargin{:});
            
            output_sizes = cellfun(@size, first_outputs, 'UniformOutput', false);
            
            % Find size of array needed to hold outputs from all windows, based on size
            % of output from a single window and on the number of windows in each
            % dimension.
            
            all_output_sizes = cellfun(@(x) [x, epoch_no], output_sizes, 'UniformOutput', false);
            
            output_dims = cellfun(@length, output_sizes);
            
            % Initialize cell holding indices for each window, and size arguments for
            % array to hold outputs from all windows.
            
            all_outputs = cellfun(@(x) nan(x), all_output_sizes, 'UniformOutput', false); % Initialize array to hold outputs from all windows.
            
            output_indices = cell(no_args_out, 1);
            
            for ao = 1:no_args_out
                
                output_indices{ao} = cell(1, output_dims(ao));
                
                output_indices{ao}(:) = {':'};
                
                all_outputs{ao}(output_indices{ao}{:}, 1) = first_outputs{ao};
                
            end
            
            for e = 1:epoch_no
                
                data = load(epoch_names{e});
                
                [outputs{1:no_args_out}] = feval(fcn_handle, data, varargin{:});
                
                for ao = 1:no_args_out
                    
                    all_outputs{ao}(output_indices{ao}{:}, e) = outputs{ao};
                    
                end
                
            end
            
            cd (present_dir)
            
            save([record_dir, sprintf('_chan%d_%s_%s', location_channels{channel_number}(s),...
                func2str(fcn_handle), make_varargin_label(varargin{:}))], 'all_outputs');
            
        end
        
    end
    
end
