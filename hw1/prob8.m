%% Problem 8

% 8. (a)
fc = 880;
fs = 1200;
I = [0:1/fs:2];
x = cos(2*pi*fc*I);
sound(x, fs);
figure;
plot(I, x);
title('signal with (fs = 1200)');
xlabel('time (s)');
ylabel('Amplitude');
saveas(gcf, 'Q8a.png');

% 8. (b)
[col, row] = size(x);
y = fft(x, row);
figure;
subplot(2,1,1);
plot(abs(y));
title('DFT amp')
xlabel('Index');
ylabel('Amplitude');
subplot(2,1,2);
plot(angle(y));
title('DFT Phase')
xlabel('Index');
ylabel('Phase');
saveas(gcf, 'Q8b.png');

N = row;
w = 2*pi * (0:(N-1)) / N;
w2 = fftshift(w);
w3 = unwrap(w2 - 2*pi);
figure;
subplot(2,1,1);
plot(w3/pi, abs(fftshift(y)));
title('DTFT Amplitude')
xlabel('radians / \pi');
ylabel('Amplitude');
subplot(2,1,2);
plot(w3/pi, angle(fftshift(y)));
title('DTFT Phase')
xlabel('radians / \pi');
ylabel('Phase');
saveas(gcf, 'Q8c.png');

y = y / fs;
figure;
subplot(2,1,1);
plot(w3/pi*fs/2, abs(fftshift(y)));
title('CTFT Amplitude')
xlabel('Freq (Hz)');
ylabel('Amplitude');
subplot(2,1,2);
plot(w3/pi*fs/2, angle(fftshift(y)));
title('CTFT Phase')
xlabel('Freq (Hz)');
ylabel('Phase');
saveas(gcf, 'Q8d.png');
