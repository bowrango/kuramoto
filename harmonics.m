% Load or generate the time-domain signal
fs = 100e3; % Hz
t = 0:1/fs:1e-2;
f0 = 9.5e3;
signal = cos(2*pi*f0*t) + 2*cos(2*pi*2*f0*t) + 3*cos(2*pi*3*f0*t);

% signal = cos(2*pi*f0*t);

N = length(signal);
Y = fft(signal);
f = (-N/2:N/2-1)*(fs/N);

% Compute magnitude in dB
magY = abs(fftshift(Y))/N; % Normalize
mag_db = 20*log10(magY);

% Plot harmonic spectrum
figure;
stem(f/fs, mag_db, 'bo'); % Normalize frequency by fs
xlabel('harmonic number');
ylabel('value db [20*log10(mag)]');
title('Harmonic Spectrum');
grid on;