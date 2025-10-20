
Fs = 1000;
L = 1e4;
t  = (0:L-1).' / Fs;

% test multi-tone signal
tones = [1;100;200;300;400];
s = sum(sin(tones*2*pi*t.')).';

% figure
% plot(t, s)
% grid on
% xlabel('Time (s)'); 
% ylabel('Amplitude');
% title('Time domain')

Fn = Fs/2;

NFFT = 2^nextpow2(L);
S = fft(s, NFFT);
% frequency axis single-sided
Fv = (0:NFFT/2) * (Fs/NFFT);

% no window: single-sided
S1 = S / L;                 
A1 = abs(S1(1:NFFT/2+1));
A1(2:end-1) = 2*A1(2:end-1); % single-sided correction

figure; 
plot(Fv, A1);
grid on;
hold on

% with Hamming window
w  = hamming(L);
cg = sum(w)/L; % coherent gain
Sw = fft(s.*w, NFFT) / (L*cg);
Aw = abs(Sw(1:NFFT/2+1));
Aw(2:end-1) = 2*Aw(2:end-1);

plot(Fv, Aw)
xlabel('Frequency (Hz)'); 
ylabel('Amplitude')
title('Amplitude Spectrum (single-sided)')
legend('No Window','Hamming Window')