Fs = 1000;
L = 1e4;
t = linspace(0, L-1, L).'/Fs;                    
s = sum(sin([1;100;200;300;400]*2*pi*t.')).';
figure
plot(t, s)
grid

Fn = Fs/2;
NFFT = 2^nextpow2(L);
FTs = fft(s,NFFT)/L;
Fv = linspace(0, 1, NFFT/2+1)*Fn;
Iv = 1:numel(Fv);

figure
plot(Fv, abs(FTs(Iv)))
grid
xlabel('Frequency')
ylabel('Magnitude')
title('Without Window')

FTs = fft(s.*hamming(L),NFFT)/L;
plot(Fv, abs(FTs(Iv)))
grid
xlabel('Frequency')
ylabel('Magnitude')
title('With Window')