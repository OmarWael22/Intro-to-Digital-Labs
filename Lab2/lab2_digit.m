

%% Reconstruction from oversampling
t=0:0.001:1;% time signal
fs = 1000; % 1/0.001
y=2*cos(2*pi*5*t);
figure;
subplot(3,1,1);
plot(t,y);
title('original signal');
%butter(order, fcut/ fs/2 , type )
[B,A] = butter(3,5/5000,'low' ); % butter fly filter
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];
% Adding zeros enhances the signal display and don't change the
%spectrum,it changes sampling freq. only
%  the length of the signal is increased by a factor of 10,
%  which effectively increases the sampling rate by a factor of 10 as well.
%  This results in a smoother waveform with higher resolution in the time domain
%  , and a higher frequency resolution in the frequency domain.
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);
subplot(3,1,2)
plot(t,filtered_signal,'r' )
xlabel('time')
title('oversampled time domian signal')
s=fft(filtered_signal);
s=fftshift(s);
fs = 10000; %1/0.001 *10 as we added 10 zeros
freq=linspace(-fs/2,fs/2,length(s));
subplot(3,1,3)
plot(freq,abs(s))
xlabel('freq')
title('freq domain oversampled signal')

%% Construction from minimum sampling

t=0:0.1:1; 
fs = 10; % fs=1/0.1
y = 2*cos(2*pi*5*t);
[B,A] = butter(10,5/50,'low' );

zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
zero_added_signal(i*10)=y(i);
end

zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);

figure
subplot(2,1,1)
plot(t,filtered_signal,'r' )
xlabel('time')
title('minimum sampled signals')

s=fft(filtered_signal);
s=fftshift(s);
fs = 100; % why 100?? 
% fs=1/0.1 * 10 as we added 10 zeros later which inreases the sampling freq by 10
freq=linspace(-fs/2,fs/2,length(s));

subplot(2,1,2)
plot(freq,abs(s))
xlabel('freq')
title('magnitude of minimum sampled signals')

%% construction from undersampling sampling

t=0:0.2:1;
fs = 5; %1/0.2
y=2*cos(2*pi*5*t);
[B,A] = butter(10,5/25,'low' );
zero_added_signal=zeros(1,length(y)*10);
for i=1:length(y)
zero_added_signal(i*10)=y(i);
end
zero_added_signal(1:9)=[];
t=linspace(0,1,length(zero_added_signal));
filtered_signal = filter(B,A,zero_added_signal);

subplot(2,1,1)
plot(t,filtered_signal,'r' )
xlabel('time')
title('undersampled signals')

s=fft(filtered_signal);
s=fftshift(s);
fs=50;
% fs=1/0.2 * 10 as we added 10 zeros later which inreases the sampling freq by 10
freq=linspace(-fs/2,fs/2,length(s));

subplot(2,1,2)
plot(freq,abs(s))
xlabel('freq')
title('magnitude of undersampled signals')


%% quantization using pi function
% Parameters
A = 1;      % Amplitude
f = 2;      % Frequency
Fs = 4000;  % Sampling frequency

% Generate sinusoidal wave
t = 0:1/Fs:1;           % Time vector
x = A*sin(2*pi*f*t);    % Sinusoidal wave

% Plotting the original sine wave to compare it to the quantized signal
figure;
plot(t,x);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Original Sine Wave');
%%
% Calculate mean square quantization error for n = 3, 4, 5, 10
n_vec = [3 4 5 10];
MSE = zeros(size(n_vec));

for i = 1:length(n_vec)
    n = n_vec(i);               % Number of bits for integer or fraction parts
    m = 2*n + 1;                % Number of bits for the whole word, including the sign bit
    xq = fi(x, 1, m, n);        % Quantize the signal using fixed-point arithmetic
    xq = double(xq);            % Convert back to double for further processing
    xq = abs(xq);
    xb = dec2bin(xq, m);         % Quantize the signal using fixed-point arithmetic
    e = x - xq;                 % Calculate quantization error
    MSE(i) = mean(e.^2);        % Calculate mean square quantization error
    
    % representation
    figure;
    plot(t,xq);
    % Add a main title to the figure
    title('Quantized signal for n = 10' );
end

figure;
plot(n_vec , MSE , 'r')
xlabel('number of bits n')
ylabel('MSE')
title('quantization using pi function')


% Note : Increasing the number of bits for quantization 
% improves the accuracy of the quantized signal,
% which results in a lower mean square quantization error.
% When n is small, 
% the number of bits used to represent the signal is limited, 
% which can result in significant quantization errors. 
% As n increases, 
% more bits are used to represent the signal, allowing for a more accurate representation of the original signal,
% and reducing the quantization error.

%% quantization using quantize function
% Parameters
A = 1;      % Amplitude
f = 2;      % Frequency
Fs = 4000;  % Sampling frequency
% Generate sinusoidal wave
t = 0:1/Fs:1;           % Time vector
x = A*sin(2*pi*f*t);    % Sinusoidal wave

% Plotting the original sine wave to compare it to the quantized signal
figure;
plot(t,x);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Original Sine Wave');
%% Calculate mean square quantization error for n = 3, 4, 5, 10
n_vec = [3 4 5 10];
MSE = zeros(size(n_vec));

for i = 1:length(n_vec)
    n = n_vec(i);
    m = 2*n + 1;                % word size
    q = quantizer([m n]);       % Create quantizer object
    xq = quantize(q, x);        % Quantize the signal
    xq = abs(xq);               % absolute of quantized signal
    xb = dec2bin(xq, m);        % Convert quantized samples to binary
    e = x - double(xq);         % Calculate quantization error
    MSE(i) = mean(e.^2);        % Calculate mean square quantization error
    
    % representation
    figure;
    plot(t,xq);
    % Add a main title to the figure
    title('Quantized signal for n = 10' );
end

figure;
plot(n_vec, MSE , 'r')  % Plot mean square quantization error
title('quantization using qunatizate function')
xlabel('Number of bits')
ylabel('Mean square quantization error')


%% non uniform quantization using compand
% Parameters
A = 1;      % Amplitude
f = 2;      % Frequency
Fs = 4000;  % Sampling frequency
% Generate sinusoidal wave
t = 0:1/Fs:1;           % Time vector
x = A*sin(2*pi*f*t);    % Sinusoidal wave
% Plotting the original sine wave to compare it to the quantized signal
figure;
plot(t,x);
xlabel('Time (s)');
ylabel('Amplitude (V)');
title('Original Sine Wave');
% Perform non-uniform quantization using compand
mu = 255;                                      % Companding law parameter
x_compand = compand(x, mu, max(x), 'mu/compressor');  % Compress signal using companding law
%% Calculate mean square quantization error for n = 3, 4, 5, 10
n_vec = [3 4 5 10];
MSE = zeros(size(n_vec));

for i = 1:length(n_vec)
    n = n_vec(i);
    m = 2*n + 1;
    q = quantizer([m n]);
    xq = quantize(q, x_compand);                           % Quantize the compressed signal
    x_expand = compand(xq, mu, max(xq),'mu/expander');           % Expand the quantized signal
    x_expand = abs(double(x_expand));     % Convert back to double and take absoltue
    xb = dec2bin(x_expand, m);    % Convert the quantized samples to binary        
    e = x - double(x_expand);             % Calculate quantization error
    MSE(i) = mean(e.^2);            % Calculate mean square quantization error
    
    % representation
    figure;
    plot(t,xq);
    % Add a main title to the figure
    title('Quantized signal for n = 10' );
end
figure
plot(n_vec, MSE , 'r')           % Plot mean square quantization error
title('quantization using compand function')
xlabel('Number of bits')
ylabel('Mean square quantization error')



