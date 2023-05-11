clear all; close all; clc;

% Parameters
nBits = 5;   % number of bits
Ts = 1e-3;       % symbol period
T = nBits*Ts;    % total time
t = 0:Ts:T-Ts;   % time vector

% Generate random bits
bits = randi([0 1],1, nBits);

% Define pulse shape
pulse_width = 10;  % Width of square pulse in samples
pulse_shape = ones(1, pulse_width);  % Square pulse

% Modulate using different line codes
nrtz = bits.*2-1;                         % non-return to zero

% Apply pulse amplitude modulation (PAM) to NRZI signal
pam_signal_nrtz = zeros(1, length(nrtz*pulse_width));
for i = 1:length(nrtz)
    start_idx = (i-1)*pulse_width+1;
    end_idx = start_idx+pulse_width-1;
    pam_signal_nrtz(start_idx:end_idx) = nrtz(i) * pulse_shape;
end

% Modulate using NRTZ inverted line code
nrzi = zeros(1, length(bits)*2);
% Set initial signal level to positive
signal_level = -1;

% Loop through binary data and encode using NRZI
j = 1; % Index variable for the AMI output vector
for i = 1:length(bits)
    if bits(i) == 0
        % No transition, keep same signal level
        nrzi(j:j+1) = signal_level;
    else
        % Transition, invert signal level
        signal_level = -signal_level;
        nrzi(j:j+1) = signal_level;
    end
     j = j + 2; % Increment the index variable by 2 for each bit
end

% Apply pulse amplitude modulation (PAM) to NRZI signal
pam_signal_nrzi = zeros(1, length(nrzi*pulse_width));
for i = 1:length(nrzi)
    start_idx = (i-1)*pulse_width+1;
    end_idx = start_idx+pulse_width-1;
    pam_signal_nrzi(start_idx:end_idx) = nrzi(i) * pulse_shape;
end

% Modulate using RTZ line code
rtz = zeros(1, length(bits)*2);
for i = 1:length(bits)
    if bits(i) == 1
        rtz(2*i-1) = 1;
    else
        rtz(2*i-1) = 0;
    end
end

% Apply pulse amplitude modulation (PAM) to NRZI signal
pam_signal_rtz = zeros(1, length(rtz*pulse_width));
for i = 1:length(rtz)
    start_idx = (i-1)*pulse_width+1;
    end_idx = start_idx+pulse_width-1;
    pam_signal_rtz(start_idx:end_idx) = rtz(i) * pulse_shape;
end

% Modulate using AMI (Alternate Mark Inversion) line code
pulse = -1;
ami = zeros(1, 2*length(bits)); % Initialize the AMI output vector with zeros
j = 1; % Index variable for the AMI output vector
for i = 1:length(bits)
    if bits(i) == 1
        if pulse == -1
            ami(j:j+1) = 1;     
        else
            ami(j:j+1) = -1;
        end
        pulse = -pulse; % Toggle the pulse
    end
    j = j + 2; % Increment the index variable by 2 for each bit
end

% Apply pulse amplitude modulation (PAM) to NRZI signal
pam_signal_ami = zeros(1, length(ami*pulse_width));
for i = 1:length(ami)
    start_idx = (i-1)*pulse_width+1;
    end_idx = start_idx+pulse_width-1;
    pam_signal_ami(start_idx:end_idx) = ami(i) * pulse_shape;
end

% Modulate using Manchester line code
MAN = zeros(1, length(bits)*2);
for i = 1:length(bits)
    if bits(i) == 1
        MAN(2*i-1) = 1;
    else
        MAN(2*i) = 1;
    end
end

% Apply pulse amplitude modulation (PAM) to NRZI signal
pam_signal_man = zeros(1, length(MAN*pulse_width));
for i = 1:length(MAN)
    start_idx = (i-1)*pulse_width+1;
    end_idx = start_idx+pulse_width-1;
    pam_signal_man(start_idx:end_idx) = MAN(i) * pulse_shape;
end

% Modulate using MLT3 line code
% define the three signal levels
levels = [-1 0 1 0];

% initialize the MLT-3 encoded signal
mlt = zeros(1, 2*length(bits));

% initialize the current signal level
current_level = 2;

% loop through each bit in the sequence
j = 1;
for i = 1:length(bits)
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
    j = j + 2 ;  % Increment the index variable by 2 for each bit
end

% Apply pulse amplitude modulation (PAM) to NRZI signal
pam_signal_mlt = zeros(1, length(mlt*pulse_width));
for i = 1:length(mlt)
    start_idx = (i-1)*pulse_width+1;
    end_idx = start_idx+pulse_width-1;
    pam_signal_mlt(start_idx:end_idx) = mlt(i) * pulse_shape;
end



%plot original binary sequence
figure;
stem(bits);xlim([-1 nBits+2]);ylim([-2 2]);title('Original Sequence');

% Plot modulated signals
figure;
subplot(6,1,1); plot(pam_signal_nrtz); axis([0 length(pam_signal_nrtz) -1.5 1.5]); ylim([-1.5 1.5]); title('Non-return to zero');
subplot(6,1,2); plot(pam_signal_nrzi); axis([0 length(pam_signal_nrzi) -1.5 1.5]);ylim([-1.5 1.5]); title('Non-return to zero inverted');
subplot(6,1,3); plot(pam_signal_rtz); axis([0 length(pam_signal_rtz) -1.5 1.5]); ylim([-1.5 1.5]); title('Return to zero');
subplot(6,1,4); plot(pam_signal_ami); axis([0 length(pam_signal_ami) -1.5 1.5]); ylim([-1.5 1.5]); title('Alternate mark inversion');
subplot(6,1,5);plot(pam_signal_man); axis([0 length(pam_signal_man) -1.5 1.5]); ylim([-1.5 1.5]); title('Manchester coding');
subplot(6,1,6); plot(pam_signal_mlt); axis([0 length(pam_signal_mlt) -1.5 1.5]); ylim([-1.5 1.5]); title('Multi-level transmission');
xlabel('Time (s)');

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


