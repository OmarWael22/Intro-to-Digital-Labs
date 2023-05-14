clear; close all; clc;

% Parameters
nBits = 5;   % number of bits
Ts = 1e-3;       % symbol period
T = nBits*Ts;    % total time
t = 0:Ts:T-Ts;   % time vector

% Generate random bits
bits = randi([0 1],1, nBits);
bits = [ 0 1 0 1 1];
%plot original binary sequence
figure;
stem(bits);xlim([-1 nBits+2]);ylim([-2 2]);title('Original Sequence');
grid on;

% Define pulse shape
pulse_width = 10;  % Width of square pulse in samples
pulse_shape = ones(1, pulse_width);  % Square pulse

% Modulate using different line codes
%%%%%% non-return to zero (NRZ)
nrtz = bits.*2-1;                         

%%%%%%% NRTI
nrzi = zeros(1, length(bits)*2);

%%%%%%% RTZ
rtz = zeros(1, length(bits)*2);
% Set initial signal level to positive
signal_level = -1;

%%%%%%% AMI
pulse = -1;
ami = zeros(1, 2*length(bits)); % Initialize the AMI output vector with zeros
j = 1; % Index variable for the AMI output vector

%%%%%%% MAN 
MAN = zeros(1, length(bits)*2);


%%%%% MLT3
% define the three signal levels
levels = [-1 0 1 0];

% initialize the MLT-3 encoded signal
mlt = zeros(1, 2*length(bits));

% initialize the current signal level of MLT3
current_level = 2;

for i = 1:length(bits)
    
    % Loop through binary data and encode using NRZI
    if bits(i) == 0
        % No transition, keep same signal level
        nrzi(j:j+1) = signal_level;
    else
        % Transition, invert signal level
        signal_level = -signal_level;
        nrzi(j:j+1) = signal_level;
    end
    
    % Modulate using RTZ line code
    if bits(i) == 1
        rtz(2*i-1) = 1;
    else
        rtz(2*i-1) = 0;
    end
    
    % Modulate using AMI (Alternate Mark Inversion) line code
    if bits(i) == 1
        if pulse == -1
            ami(j:j+1) = 1;     
        else
            ami(j:j+1) = -1;
        end
        pulse = -pulse; % Toggle the pulse
    end
    
    % Modulate using Manchester line code
    if bits(i) == 1
        MAN(2*i-1) = 1;
    else
        MAN(2*i) = 1;
    end
    
    % Modulate using MLT3 line code
    % if the current bit is 0, the signal level remains the same
    if bits(i) == 0
        mlt(j:j+1) = levels(current_level);
    % if the current bit is 1, the signal level changes to the next available level
    else
        % change the signal level
        current_level = mod(current_level, 4) + 1;
        % set the new signal level for the current bit
        mlt(j:j+1) = levels(current_level);
    end
     j = j + 2; % Increment the index variable by 2 for each bit
end

% Signals to PAM signals To make plotting like Square waves
signals = {nrtz, nrzi, rtz, ami, MAN, mlt};
signals_names = {'NRTZ', 'NRTZI', 'RTZ', 'AMI', 'MAN', 'MLT3'};
pam_signals = cell(size(signals));

for k = 1:numel(signals)
    signal = signals{k};
    pam_signal = zeros(1, length(signal*pulse_width));
    for i = 1:length(signal)
        start_idx = (i-1)*pulse_width+1;
        end_idx = start_idx+pulse_width-1;
        pam_signal(start_idx:end_idx) = signal(i) * pulse_shape;
    end
    pam_signals{k} = pam_signal;
end

% plotting signals
figure;
for k = 1:numel(pam_signals)
    subplot(numel(pam_signals), 1, k);
    plot(pam_signals{k});
    axis([0 length(pam_signals{k}) -1.5 1.5]); ylim([-1.5 1.5]);
    title(sprintf('PAM Signal for %s', upper(signals_names{k})));
    xlabel('Time (s)');
end

% Calculate power spectral density
[P_nrtz, f_nrtz] = periodogram(nrtz,[],[],1/Ts,'centered');
[P_nrzi, f_nrzi] = periodogram(nrzi,[],[],1/Ts,'centered');
[P_rtz, f_rtz] = periodogram(rtz,[],[],1/Ts,'centered');
[P_ami, f_ami] = periodogram(ami,[],[],1/Ts,'centered');
[P_manchester, f_manchester] = periodogram(MAN,[],[],1/Ts,'centered');
[P_mlt, f_mlt] = periodogram(mlt,[],[],1/Ts,'centered');

% Plot power spectral density
figure;
subplot(6,1,1); plot(f_nrtz,10*log10(P_nrtz)); title('Non-return to zero');
subplot(6,1,2); plot(f_nrzi,10*log10(P_nrzi)); title('Non-return to zero inverted');
subplot(6,1,3); plot(f_rtz,10*log10(P_rtz)); title('Return to zero');
subplot(6,1,4); plot(f_ami,10*log10(P_ami)); title('Alternate mark inversion');
subplot(6,1,5); plot(f_manchester,10*log10(P_manchester)); title('Manchester coding');
subplot(6,1,6); plot(f_mlt,10*log10(P_mlt)); title('Multi-level transmission');
xlabel('Frequency (Hz)');
ylabel('Power spectral density (dB/Hz)');

