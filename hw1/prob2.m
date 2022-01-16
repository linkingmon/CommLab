%% Problem 2

% Read audio files
[x, fs] = audioread('handel.ogg');

% 2. (a)
xt = circshift(x, 100000);
figure;
manip(x, fs, 1, 'Original Waveform');
manip(xt, fs, 2, 'Waveform with Cirshift 100000');
saveas(gcf, 'Q2a.png');

% 2. (b)
xt = flipud(x);
figure;
manip(x, fs, 1, 'Original Waveform');
manip(xt, fs, 2, 'Waveform with Time reversal');
saveas(gcf, 'Q2b.png');

% 2. (c)
xt = [x;xt];
figure;
manip(x, fs, 1, 'Original Waveform');
manip(xt, fs, 2, 'Waveform with forward-then-backward');
saveas(gcf, 'Q2c.png');

% 2. (d)
xt = upsample(x, 3);
figure;
manip(x, fs, 1, 'Original Waveform');
manip(xt, fs, 2, 'Waveform with upsample 3');
saveas(gcf, 'Q2d1.png');

xt = downsample(x, 6);
figure;
manip(x, fs, 1, 'Original Waveform');
manip(xt, fs, 2, 'Waveform with downsample 3');
saveas(gcf, 'Q2d2.png');

function manip(x, fs, num, tit)
    [row, col] = size(x);
    time = [0:1/fs:(row-1)/fs]';
    subplot(2,1,num);
    plot(time, x);
    title(tit);
    xlabel('time (s)');
    ylabel('Amplitude');
end