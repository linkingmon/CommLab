%% problem 1

% (a)
fs = 400;
I = [0:1/fs:1-1/fs];
x1 = cos(2*pi*10*I);
x2 = cos(2*pi*25*I);
x3 = cos(2*pi*50*I);
x4 = cos(2*pi*100*I);
x5 = I * 0;
x = [x1 x2 x3 x4 x5];
It = [0:1/fs:5-1/fs];


O = ones(1,fs);
Z = zeros(1,fs);
w = [O Z Z Z Z];
CTFT_draw(x.*w, fs, It);
saveas(gcf, 'Q1a.png');

% (b)
CTFT_draw(x.*circshift(w,5*fs/2), fs, It);
saveas(gcf, 'Q1b.png');

% (c)
% STFT() freq start from negative, while spectrogram() start from 0
figure;
w = ones(1,400);
noverlap = 200;
Nfft = 512;
[s1, f, t] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t, f, flipud(abs(s1)));
set(gca,'YDir','normal');
title('my STFT (400, 200)')
colorbar;
saveas(gcf, 'Q1d_1.png');

figure;
w = ones(1,400);
noverlap = 200;
Nfft = 512;
[s2, f2, t2] = spectrogram(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('matlab spectrogram (400, 200)')
colorbar;
saveas(gcf, 'Q1d_2.png');

% (d) adjust window size, noverlap
figure;
w = ones(1,100);
noverlap = 50;
Nfft = 512;
[s1, f, t] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t, f, flipud(abs(s1)));
set(gca,'YDir','normal');
title('my STFT (100, 50)')
colorbar;
saveas(gcf, 'Q1d_3.png');

figure;
w = ones(1,100);
noverlap = 50;
Nfft = 512;
[s2, f2, t2] = spectrogram(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('matlab spectrogram (100, 50)')
colorbar;
saveas(gcf, 'Q1d_4.png');

figure;
w = ones(1,400);
noverlap = 50;
Nfft = 512;
[s1, f, t] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t, f, flipud(abs(s1)));
set(gca,'YDir','normal');
title('my STFT (400, 50)')
colorbar;
saveas(gcf, 'Q1d_5.png');

figure;
w = ones(1,400);
noverlap = 50;
Nfft = 512;
[s2, f2, t2] = spectrogram(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('matlab spectrogram (400, 50)')
colorbar;
saveas(gcf, 'Q1d_6.png');

% (e)
[x, fs] = audioread('handel.ogg');
figure;
w = ones(1,100);
noverlap = 50;
Nfft = 128;
[s2, f2, t2] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('original spectrogram')
colorbar;
saveas(gcf, 'Q1e_1.png');


[x, fs] = audioread('handel_square.wav');
figure;
w = ones(1,100);
noverlap = 50;
Nfft = 128;
[s2, f2, t2] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('square spectrogram')
colorbar;
saveas(gcf, 'Q1e_2.png');


[x, fs] = audioread('handel_4-bit_Quantization.wav');
figure;
w = ones(1,100);
noverlap = 50;
Nfft = 128;
[s2, f2, t2] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('quantization spectrogram')
colorbar;
saveas(gcf, 'Q1e_3.png');


[x, fs] = audioread('handel_hardlimit_0.1.wav');
figure;
w = ones(1,100);
noverlap = 50;
Nfft = 128;
[s2, f2, t2] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('hardlimit spectrogram')
colorbar;
saveas(gcf, 'Q1e_4.png');

[x, fs] = audioread('handel_modulation_1000.wav');
figure;
w = ones(1,100);
noverlap = 50;
Nfft = 128;
[s2, f2, t2] = STFT(x, w, noverlap, Nfft, fs);
imagesc(t2, f2, abs(s2)/fs);
set(gca,'YDir','normal');
title('modulation spectrogram')
colorbar;
saveas(gcf, 'Q1e_5.png');
clear all;

function [W, A] = CTFT_draw(x, fs, It)   
    figure;
    subplot(3,1,1);
    plot(It, x);
    ylabel('Amp');
    xlabel('Time');
    
    [col, row] = size(x);
    y = fft(x, row);

    N = row;
    w = 2*pi*(0:(N-1)) / N;
    w2 = fftshift(w);
    w3 = unwrap(w2 - 2*pi);

    y = y / fs;
    W = w3/pi*fs/2;
    A = abs(fftshift(y));
    subplot(3,1,2);
    plot(w3/pi*fs/2, abs(fftshift(y)));
    title('CTFT Amplitude');
    xlabel('Freq (Hz)');
    ylabel('Amplitude');
    subplot(3,1,3);
    plot(w3/pi*fs/2, angle(fftshift(y)));
    title('CTFT Phase');
    xlabel('Freq (Hz)');
    ylabel('Phase');
end

function [W, A] = CTFT(x, fs, Nfft)   
    y = fft(x, Nfft);
    w = 2*pi*(0:(Nfft-1)) / Nfft;
    w2 = fftshift(w);
    w3 = unwrap(w2 - 2*pi);

    y = y / fs;
    W = w3/pi*fs/2;
    A = fftshift(y);
end

function [S, f, t] = STFT(x, window, Noverlap, Nfft, fs)  
    % If x is a signal of length Nx, then s has k columns, where
    % k = ⌊(Nx – noverlap)/(length(window) – noverlap)⌋ if window is a vector.
    % If x is real and nfft is even, then s has (nfft/2 + 1) rows.
    % If x is real and nfft is odd, then s has (nfft + 1)/2 rows.
    % If x is complex, then s has nfft rows.
    assert(length(window) > Noverlap, 'size of window must be larger than Noverlap')
    col_size = Nfft;
    row_size = floor( (length(x) - Noverlap) / (length(window) - Noverlap) );
    S = zeros(col_size, row_size);
    t = zeros(1, row_size);
    from = 1;
    to = length(window);
    for i = 1 : row_size
        X = x(from:to);
        [f, A] = CTFT(X, fs, Nfft);
        S(1:end,i) = A(1:col_size);
        t(1, i) = ( from + floor(length(window) / 2) - 1)/ fs;
        from = from + (length(window) - Noverlap);
        to = to + (length(window) - Noverlap);
    end
    f = f(1:col_size);
end
