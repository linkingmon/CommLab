%% Problem 1

% Read audio files
[x, fs] = audioread('handel.ogg', 'native');
whos x
sound(x, fs);
% changing sampling freq and plot
figure;
manip(x, fs, 1, 'Original Waveform fs = 44100');
fst = 20000;
manip(x, fst, 2, 'Waveform with fs = 20000');
fst = 80000;
manip(x, fst, 3, 'Waveform with fs = 80000');
saveas(gcf, 'Q1.png');

function manip(x, fs, num, tit)
    [row, col] = size(x);
    time = [0:1/fs:(row-1)/fs]';
    subplot(3,1,num);
    plot(time, x);
    title(tit);
    xlabel('time (s)');
    ylabel('Amplitude');
end