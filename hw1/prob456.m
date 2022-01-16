%% Problem 4 5 6

% Read audio files
[x, fs] = audioread('handel.ogg');

% 4. (a) (b) (c) (d) change bit
% assume amp is in [-1, 1]
xmax = 1;
bit = 4;
bit2 = 6;
level = 2^bit;
level2 = 2^bit2;
xt = quantizer_L_level(x, xmax, level)';
xt2 = quantizer_L_level(x, xmax, level2)';
% sound(xt2, fs);
figure;
manip(xt, fs, 1, 'Waveform with 4-bit Quantization');
audiowrite('handel_4-bit_Quantization.wav', xt, fs)
manip(xt2, fs, 2, 'Waveform with 6-bit Quantization');
saveas(gcf, 'Q4.png');

% 5. 
fc = 1000;
I = [1:size(x)];
clear i;
xt = x.*exp(i*2*pi*fc*I'/fs);
% sound(real(xt), fs);
figure;
manip2(abs(xt), fs, 1, 'Amp with modulation fc = 100', 'Amp');
audiowrite('handel_modulation_100.wav', abs(xt), fs)
manip2(angle(xt), fs, 2, 'Phase with modulation fc = 100', 'Phase');
saveas(gcf, 'Q5.png');

% 6. (Noise)
var = 0.01;
xt = x + sqrt(var)*randn(size(x));
sound(xt, fs);
figure;
manip(x, fs, 1, 'Original waveform');
manip(xt, fs, 2, 'Waveform with Gaussian noise with var = 0.01');
saveas(gcf, 'Q6.png');

function y = quantizer_L_level(x, xmax, level)
    delta = 2 * xmax / level;
    partition = [-xmax:delta:xmax];
    codebook = [0,-(level-1)*delta/2:delta:(level-1)*delta/2,0];
    [I, y] = quantiz(x,partition,codebook); 
end


function manip(x, fs, num, tit)
    [row, col] = size(x);
    time = [0:1/fs:(row-1)/fs]';
    subplot(2,1,num);
    plot(time, x);
    title(tit);
    xlabel('time (s)');
    ylabel('Amplitude');
end

function manip2(x, fs, num, tit, yl)
    [row, col] = size(x);
    time = [0:1/fs:(row-1)/fs]';
    subplot(2,1,num);
    plot(time, x);
    title(tit);
    xlabel('time (s)');
    ylabel(yl);
end