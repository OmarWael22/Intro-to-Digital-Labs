clear all; close all; clc;

% Parameters
nBits = 5;   % number of bits
Ts = 1e-3;       % symbol period
T = nBits*Ts;    % total time
t = 0:Ts:T-Ts;   % time vector

% Generate random bits
bits = randi([0 1],1, nBits);

% Modulate using different line codes
nrtz = bits.*2-1;                         % non-return to zero
nrzi = 1-2*bits;                          % non-return to zero inverted

% Modulate using RTZ line code
rtz = zeros(1, length(bits)*2);
for i = 1:length(bits)
    if bits(i) == 1
        rtz(2*i-1) = 1;
    else
        rtz(2*i-1) = 0;
    end
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

% Modulate using Manchester line code
MAN = zeros(1, length(bits)*2);
for i = 1:length(bits)
    if bits(i) == 1
        MAN(2*i-1) = 1;
    else
        MAN(2*i) = 1;
    end
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

%plot original binary sequence
figure;
stem(bits);xlim([-1 nBits+2]);ylim([-2 2]);title('Original Sequence');

% Plot modulated signals
figure;
subplot(6,1,1); plot(nrtz); ylim([-1.5 1.5]); title('Non-return to zero');
subplot(6,1,2); plot(nrzi); ylim([-1.5 1.5]); title('Non-return to zero inverted');
subplot(6,1,3); plot(rtz); ylim([-1.5 1.5]); title('Return to zero');
subplot(6,1,4); plot(ami); ylim([-1.5 1.5]); title('Alternate mark inversion');
subplot(6,1,5); plot(MAN); ylim([-1.5 1.5]); title('Manchester coding');
subplot(6,1,6); plot(mlt); ylim([-1.5 1.5]); title('Multi-level transmission');
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


