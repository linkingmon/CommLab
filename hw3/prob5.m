%% problem 2

% assume amp is in [-1, 1]
clear all;
[x, fs] = audioread('handel.ogg');
xmax = 1;
bit = 4;
level = 2^bit;
xt = quantizer_L_level(x, xmax, level)';
y = pcm_enc(xt, 4);
% length_y = length(y);
num1 = 16;
num2 = 16;
num3 = 16;
constellation1 = 'PAM';
constellation2 = 'PSK';
constellation3 = 'QAM';
mapping1 = 'gray';
mapping2 = 'gray';
mapping3 = 'gray';
y_rec2 = Run(y, num1, 1, constellation1, mapping1, 'MD', 20);
y_rec3 = Run(y, num1, 1, constellation1, mapping1, 'MD', 10);
y_rec4 = Run(y, num1, 1, constellation1, mapping1, 'MD', 0);
y_rec6 = Run(y, num2, 1, constellation2, mapping2, 'MD', 20);
y_rec7 = Run(y, num2, 1, constellation2, mapping2, 'MD', 10);
y_rec8 = Run(y, num2, 1, constellation2, mapping2, 'MD', 0);
y_rec10 = Run(y, num3, 1, constellation3, mapping3, 'MD', 20);
y_rec11 = Run(y, num3, 1, constellation3, mapping3, 'MD', 10);
y_rec12 = Run(y, num3, 1, constellation3, mapping3, 'MD', 0);

x_rec2 = pcm_dec(y_rec2, 4);
x_rec3 = pcm_dec(y_rec3, 4);
x_rec4 = pcm_dec(y_rec4, 4);
x_rec6 = pcm_dec(y_rec6, 4);
x_rec7 = pcm_dec(y_rec7, 4);
x_rec8 = pcm_dec(y_rec8, 4);
x_rec10 = pcm_dec(y_rec10, 4);
x_rec11 = pcm_dec(y_rec11, 4);
x_rec12 = pcm_dec(y_rec12, 4);

% delta = 2 * xmax / level;
% symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
% p = histc(xt, symbols);
% p = p / sum(p);
% dict = huffmandict(symbols, p);
% y = huffmanenco(xt, dict);
% length_y = length(y)
% y_rec = Run(y, 16, 1, 'QAM', 'Gray', 'MD', 10);
% xd = huffmandeco(y_rec, dict);

sound(x_rec1, fs);

% clear all;

function y = quantizer_L_level(x, xmax, level)
    delta = 2 * xmax / level;
    partition = [-xmax:delta:xmax];
    codebook = [0,-(level-1)*delta/2:delta:(level-1)*delta/2,0];
    [I, y] = quantiz(x,partition,codebook); 
end


function binary_sequence_rec = Run(binary_sequence, M, d, constellation, mapping, decision_rule, SNR)
    tit = strcat(num2str(M), '-', constellation, '-', mapping);
    A = {};
    load(strcat('cell_prob1/', tit,'.mat'), 'A');
    B = [A{:}];
    Es = sum(abs(B).^2)/length(B);
    noise_var = Es/(10^(SNR/10));
    pd = makedist('Normal');
    symbol_sequence = symbol_mapper(binary_sequence, M, d, constellation, mapping);
    save('Noise.mat','noise_var');
    NI = random(pd,1,length(symbol_sequence))*sqrt(noise_var/2);
    NQ = random(pd,1,length(symbol_sequence))*sqrt(noise_var/2);
    symbol_sequence_rec = symbol_sequence + NI + i*NQ;
    binary_sequence_rec = symbol_demapper(symbol_sequence_rec, M, d, constellation, mapping, decision_rule);
    binary_seq_len = length(binary_sequence);
    berr = sum(binary_sequence ~= binary_sequence_rec) / binary_seq_len;
    fprintf("%13s Bit error rate   %f\n", tit, berr);
    tit = strcat(tit, '-', decision_rule);
end

function binary_sequence = symbol_demapper(symbol_sequence, M, d, constellation, mapping, decision_rule)
    switch decision_rule
        case 'MD'
            C = cell(length(symbol_sequence),1);
            tit = strcat(num2str(M), '-', constellation, '-', mapping);
            A = {};
            load(strcat('cell_prob1/', tit,'.mat'), 'A');
            B = [A{:}];
            sym_len = length(symbol_sequence);
            for ii = 1:sym_len
                [minvalue, index_of_min] = min(abs(symbol_sequence(ii) - B));
                C{ii,1} = index_of_min-1;
            end
            binary_sequence = reshape(fliplr(de2bi([C{:}],log2(M)))',1,[]);
        case 'ML'
            noise_var = 0;
            load('Noise.mat', 'noise_var');
            C = cell(length(symbol_sequence),1);
            tit = strcat(num2str(M), '-', constellation, '-', mapping);
            A = {};
            load(strcat('cell_prob1/', tit,'.mat'), 'A');
            B = [A{:}];
            P = ones(1,length(B))/length(B);
            sym_len = length(symbol_sequence);
            for ii = 1:sym_len
                [maxvalue, index_of_max] = max(exp(-abs(B - symbol_sequence(ii)).^2/noise_var).*P);
                C{ii,1} = index_of_max-1;
            end
            binary_sequence = reshape(fliplr(de2bi([C{:}],log2(M)))',1,[]);
        case 'MAP'
            P = 0;
            noise_var = 0;
            load('Prob.mat','P');
            load('Noise.mat', 'noise_var');
            C = cell(length(symbol_sequence),1);
            tit = strcat(num2str(M), '-', constellation, '-', mapping);
            A = {};
            load(strcat('cell_prob1/', tit,'.mat'), 'A');
            B = [A{:}];
            sym_len = length(symbol_sequence);
            for ii = 1:sym_len
                [maxvalue, index_of_max] = max(exp(-abs(B - symbol_sequence(ii)).^2/noise_var).*P);
                C{ii,1} = index_of_max-1;
            end
            binary_sequence = reshape(fliplr(de2bi([C{:}],log2(M)))',1,[]);
        otherwise
            error("Only 'MD' 'ML' 'MAP' are available")
    end
end

function symbol_sequence = symbol_mapper(binary_sequence, M, d, constellation, mapping);
    seq = reshape(binary_sequence,log2(M),[])';
    C = cell(length(seq),1);
    tit = strcat(num2str(M), '-', constellation, '-', mapping);
    A = {};
    load(strcat('cell_prob1/', tit,'.mat'), 'A');
    B = [A{:}];
    for ii = 1:length(seq)
        C{ii,1} = B(bi2de(fliplr(seq(ii,:)))+1);
    end
    symbol_sequence = [C{:}];
    P = zeros(1,length(B));
    for ii = 1:length(B)
        P(ii) = sum(symbol_sequence(:) == B(ii));
    end
    P = P / length(symbol_sequence);
    save('Prob.mat', 'P');
end

function y = pcm_enc(x, numBits)
    level = 2^numBits;
    delta = 2 / level;
    x = (x + (level-1)*delta/2) / delta + 1;
    y = reshape(fliplr(de2bi(x))',1,[]);
end

function x = pcm_dec(y, numBits)
    level = 2^numBits;
    delta = 2 / level;
    x = bi2de(fliplr(reshape(y, 4, [])'));
    x = (x - 1)*delta - (level-1)*delta/2;
end