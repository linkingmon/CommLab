%% Problem 3

% Read audio files
[x, fs] = audioread('handel.ogg');

% 3. (a)
y = x > 0.1 | x < -0.1;
xt = x.* (1-y) + 0.1*y;
figure;
manip(xt, fs, 1, 'Waveform with Hard-limit 0.1');
audiowrite('handel_hardlimit_0.1.wav', xt, fs)
xt = x.^2;
manip(xt, fs, 2, 'Waveform with square');
audiowrite('handel_square.wav', xt, fs)
xt = -x;
manip(xt, fs, 3, 'Waveform with negation');

% 3. (b)
thres = 0.25;
y = x > thres | x < -thres;
xt = x.* (1-y) + thres*y;
manip(xt, fs, 4, 'Waveform with Hard-limit 0.25');
saveas(gcf, 'Q3.png');


function manip(x, fs, num, tit)
    [row, col] = size(x);
    time = [0:1/fs:(row-1)/fs]';
    subplot(2,2,num);
    plot(time, x);
    title(tit);
    xlabel('time (s)');
    ylabel('Amplitude');
end