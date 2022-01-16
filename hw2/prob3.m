%% problem 3(f)

% assume amp is in [-1, 1]
[x, fs] = audioread('handel.ogg');
[x1, fs] = audioread('handel_hardlimit_0.1.wav');
[x2, fs] = audioread('handel_modulation_100.wav');
[x3, fs] = audioread('handel_square.wav');
xmax = max(abs(x))+1e-4;
xmax1 = max(abs(x1))+1e-4;
xmax2 = max(abs(x2))+1e-4;
xmax3 = max(abs(x3))+1e-4;

bit = 4;
level = 2^bit;
xt = quantizer_L_level(x, xmax, level)';
xt1 = quantizer_L_level(x1, xmax1, level)';
xt2 = quantizer_L_level(x2, xmax2, level)';
xt3 = quantizer_L_level(x3, xmax3, level)';

% huffman dict
delta = 2 * xmax / level;
symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
p = histc(xt, symbols);
p = p / sum(p);
dict = huffmandict(symbols, p);


delta = 2 * xmax1 / level;
symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
p1 = histc(xt1, symbols);
p1 = p1 / sum(p1);
dict1 = huffmandict(symbols, p1);


delta = 2 * xmax2 / level;
symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
p2 = histc(xt2, symbols);
p2 = p2 / sum(p2);
dict2 = huffmandict(symbols, p2);


delta = 2 * xmax3 / level;
symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
p3 = histc(xt3, symbols);
p3 = p3 / sum(p3);
dict3 = huffmandict(symbols, p3);

% huffman encode
% original
y = huffmanenco(xt, dict);
length_y = length(y)
% hardlimit
y1 = huffmanenco(xt1, dict1);
length_y1 = length(y1)
% modulation
y2 = huffmanenco(xt2, dict2);
length_y2 = length(y2)
% square
y3 = huffmanenco(xt3, dict3);
length_y3 = length(y3)
% compactness
figure;
subplot(2,2,1);
histogram(abs(x));
title('Original')
subplot(2,2,2);
histogram(abs(x1));
title('Hardlimit')
subplot(2,2,3);
histogram(abs(x2));
title('Modulation')
subplot(2,2,4);
histogram(abs(x3));
title('Square')

% clear all;

function y = quantizer_L_level(x, xmax, level)
    delta = 2 * xmax / level;
    partition = [-xmax:delta:xmax];
    codebook = [0,-(level-1)*delta/2:delta:(level-1)*delta/2,0];
    [I, y] = quantiz(x,partition,codebook); 
end
