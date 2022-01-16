%% Problem 7

% Read audio files
[x, fs] = audioread('handel.ogg');

% 7. (Filtering)
W = 50;
lower_freq = 94 / fs;
higher_freq = 142 / fs;
h = fir1(W, [lower_freq, higher_freq]);
impz(h);
saveas(gcf, 'Q7a.png');
% freqz(h,1,512)
xt = filter(h, 1, x);
sound(xt, fs);
figure;
manip(x, fs, 1, 'Original waveform');
manip(xt, fs, 2, 'Waveform after filtering');
saveas(gcf, 'Q7b.png');


function manip(x, fs, num, tit)
    [row, col] = size(x);
    time = [0:1/fs:(row-1)/fs]';
    subplot(2,1,num);
    plot(time, x);
    title(tit);
    xlabel('time (s)');
    ylabel('Amplitude');
end


