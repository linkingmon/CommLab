%% problem 4

% 1 a 0.3
% 2 b 0.2
% 3 c 0.2
% 4 d 0.1
% 5 e 0.05
% 6 f 0.05
% 7 g 0.05
% 8 h 0.05
n = 100;
symbols = [1:8];
prob = [0.3 0.2 0.2 0.1 0.05 0.05 0.05 0.05];
dict = huffmandict(symbols, prob);
len = [];
for i = 1 : length(dict)
    len = [len length(dict{i,2}) ];
end

entropy = -sum(prob.*log2(prob))

avg_length = sum(prob.*len)

v = randsrc(n, 1, [symbols; prob]);
v_enc = huffmanenco(v, dict);
v_enc_length = length(v_enc)

sum = 0;
R = 1000;
for i = 1 : R
    v = randsrc(n, 1, [symbols; prob]);
    v_enc = huffmanenco(v, dict);
    sum = sum + length(v_enc);
end
avg_length = sum / R

avg_length_d100 = avg_length / n

n = 1000;
R = 10000;
for i = 1 : R
    v = randsrc(n, 1, [symbols; prob]);
    v_enc = huffmanenco(v, dict);
    sum = sum + length(v_enc);
end
avg_length = sum / R

avg_length_d1000 = avg_length / n

clear all;