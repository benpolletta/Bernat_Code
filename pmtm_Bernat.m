function [spec, BP] = pmtm_Bernat(data)

sampling_freq=1000;
seconds_per_epoch=4096/250;
signal_length=sampling_freq*seconds_per_epoch;

f=sampling_freq*[0:signal_length/2]/(signal_length);

data=detrend(data);

spec=pmtm(data);

line_noise_limits=[59 61; 118 122];

[no_stops,~]=size(line_noise_limits);

stop_indices=cell(no_stops,1);
for i=1:no_stops
    stop_indices{i}=find(f>=line_noise_limits(i,1) & f<=line_noise_limits(i,2));
end

for i=1:no_stops
    
    spec(stop_indices{i})=nan;
    
end

band_limits = [1 5; 5 10; 20 40; 40 90; 120 160; 160 200];
% band_labels = {'delta', 'theta', 'low-gamma', 'high-gamma', 'HFOs', 'spindle'};

% band_limits=[20 50; 50 90; 90 120];
% band_labels={'low-gamma','middle-gamma','high-gamma'};
% band_labels_long={'Low Gamma','Middle Gamma','High Gamma'};

% band_limits=[.1 4; 4 8; 10 13; 13 20; 20 50; 50 90; 90 120; 125 175];
% band_labels={'delta','theta','alpha','low-beta','beta-low-gamma','mid-gamma','high-gamma','HFOs'};

[no_bands,~]=size(band_limits);

band_indices=cell(no_bands,1);
% band_freq_labels=cell(no_bands,1);
for i=1:no_bands
    band_indices{i}=find(f>=band_limits(i,1) & f<=band_limits(i,2));
    % band_freq_labels{i}=[band_labels{i},num2str(band_limits(i,1)),'-',num2str(band_limits(i,2))];
end

for i=1:no_bands
    
    BP(i)=nansum(spec(band_indices{i}));
    
end