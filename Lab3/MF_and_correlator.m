clc ; clear ; close all;

% Simulation Parameters
num_bits = 1e5;             % number of bits
snr_range = 0:2:30;         % SNR range in dB
num_samples = 20;           % samples per waveform
sampling_instant = 20;      % sampling instant
receiver_type = "correlator";  % receiver type: matched filter or correlator
s1_amp = 1;                 % amplitude of rectangular signal s1(t)
s2_amp = 0;                 % amplitude of zero signal s2(t)
 s1_waveform=ones(1,20);
 s2_waveform=zeros(1,20);

% Generate random binary data
bits = randi([0 1],1,num_bits);  

% Represent each bit with waveform
 waveform = reshape(repmat(bits,20,1),1,[]);


% Apply noise to waveform
for snr_idx = 1:length(snr_range)
    snr = snr_range(snr_idx);
    snr_lin = 10^(snr/10);  % convert dB to linear scale
    
    % calculate noise power based on SNR
    s1_power = s1_amp^2;         
    noise_power = s1_power/snr_lin;
    
    % add noise to waveform
    noisy_waveform = awgn(waveform,snr,'measured','linear');  
    
    % Apply convolution process in receiver
   
    
    if receiver_type == "matched"
           % matched filter
            filter=flip(s1_waveform-s2_waveform);
           detected_bits=bits;
           for i=0:length(bits)-1
               frame=( noisy_waveform(i*20+1:(i+1)*20) );
            output_samples = conv(frame,filter);
            detected_bits(i+1) = output_samples(round((length(frame)+length(filter)-1)/2 ));
           end
        
    elseif receiver_type == "correlator"
        filter = ones(1,length(noisy_waveform));    % correlator
        output_samples = filter.*noisy_waveform;
        detected_bits = sum(reshape(output_samples,num_samples,[]));
        
    else
        error("Invalid receiver type");
    end
    
    % Decide whether the Rx_sequence is ‘1’ or ‘0’ by comparing with threshold
    threshold = mean(detected_bits);
    detected_bits = detected_bits > threshold;
    
    % Calculate bit error rate (BER)
    num_errors = nnz(xor(bits,detected_bits));
    ber(snr_idx) = num_errors/num_bits;
end

% Plot BER vs SNR curve
semilogy(snr_range,ber,'-o');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR Curve');
grid on;
