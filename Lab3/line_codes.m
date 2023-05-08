clear all; close all; clc;

% Parameters
nBits = 10000;   % number of bits
Ts = 1e-3;       % symbol period
T = nBits*Ts;    % total time
t = 0:Ts:T-Ts;   % time vector

% Generate random bits
bits = randi([0 1],1,nBits);

% Modulate using different line codes
nrtz = bits.*2-1;                         % non-return to zero
nrzi = 1-2*bits;                          % non-return to zero inverted
rtz = [bits; -bits]; rtz = rtz(:).';      % return to zero
ami = 2*bits-1; ami(ami==0) = -1;         % alternate mark inversion
manchester = 2*bits-1; manchester(2:2:end) = -manchester(2:2:end); % manchester
mlt = zeros(1,length(bits)*2);            % multi-level transmission
mlt(1:2:end) = bits*2-1;
mlt(2:2:end) = -bits*2+3;

% Plot modulated signals
figure;
subplot(6,1,1); plot(nrtz); ylim([-1.5 1.5]); title('Non-return to zero');
subplot(6,1,2); plot(nrzi); ylim([-1.5 1.5]); title('Non-return to zero inverted');
subplot(6,1,3); plot(rtz); ylim([-1.5 1.5]); title('Return to zero');
subplot(6,1,4); plot(ami); ylim([-1.5 1.5]); title('Alternate mark inversion');
subplot(6,1,5); plot(manchester); ylim([-1.5 1.5]); title('Manchester coding');
subplot(6,1,6); plot(mlt); ylim([-1.5 1.5]); title('Multi-level transmission');
xlabel('Time (s)');

% Calculate power spectral density
[P_nrtz, f_nrtz] = periodogram(nrtz,[],[],1/Ts,'centered');
[P_nrzi, f_nrzi] = periodogram(nrzi,[],[],1/Ts,'centered');
[P_rtz, f_rtz] = periodogram(rtz,[],[],1/Ts,'centered');
[P_ami, f_ami] = periodogram(ami,[],[],1/Ts,'centered');
[P_manchester, f_manchester] = periodogram(manchester,[],[],1/Ts,'centered');
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
