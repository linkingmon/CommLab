%% problem 2

% assume amp is in [-1, 1]
clear all;
[x, fs] = audioread('handel.ogg');
xmax = 1;
bit = 6;
level = 2^bit;

% quantizer
fprintf("Quantizing...\n");
xt = quantizer_L_level(x, xmax, level)';

% huffman encode
fprintf("Huffman encoding...\n");
delta = 2 * xmax / level;
symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
p = histc(xt, symbols);
p = p / sum(p);
dict = huffmandict(symbols, p);
y_huffen = huffmanenco(xt, dict);

% convolution encode
fprintf("ECC encoding...\n");
impulse_response = [1 1 1 1; 1 1 0 1; 0 1 0 1];
y_conven = convolutional_enc(y_huffen, impulse_response);

% constellation mapping & demapping
fprintf("Passing constellation...\n");
num1 = 8;
constellation1 = 'PAM';
mapping1 = 'gray';
SNR = 5;
y_cons = Run(y_conven, num1, 1, constellation1, mapping1, 'MD', SNR);

% convolution decode
fprintf("ECC decoding...\n");
y_convde = convolutional_dec(y_cons, impulse_response);

% huffman decode
fprintf("Huffman decoding...\n");
y_huffde = huffmandeco(y_convde, dict);

function y = quantizer_L_level(x, xmax, level)
    delta = 2 * xmax / level;
    partition = [-xmax:delta:xmax];
    codebook = [0,-(level-1)*delta/2:delta:(level-1)*delta/2,0];
    [I, y] = quantiz(x,partition,codebook); 
end

function binary_sequence_rec = Run(binary_sequence, M, d, constellation, mapping, decision_rule, SNR)
    tit = strcat(num2str(M), '-', constellation, '-', mapping);
    A = {};
    load(strcat('../hw3/cell_prob1/', tit,'.mat'), 'A');
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
            load(strcat('../hw3/cell_prob1/', tit,'.mat'), 'A');
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
            load(strcat('../hw3/cell_prob1/', tit,'.mat'), 'A');
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
            load(strcat('../hw3/cell_prob1/', tit,'.mat'), 'A');
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
    load(strcat('../hw3/cell_prob1/', tit,'.mat'), 'A');
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

function encoded_data = convolutional_enc(binary_data, impulse_response)
    [row, col] = size(impulse_response);
    len = length(binary_data);
    for ii = 1 : row
        encoded_data(ii,:) = mod(conv(binary_data, impulse_response(ii,:)), 2);
    end
    encoded_data = reshape(encoded_data,1,[]);
end

function decoded_data = convolutional_dec(binary_data, impulse_response)
    [row, col] = size(impulse_response);
    num_state = 2^(col-1);
    states = de2bi([0:num_state-1], col-1, 'left-msb');
    Cu = zeros(num_state, row);
    Cd = zeros(num_state, row);
    for ii = 1 : row
        Cu(:,ii) = mod(sum([states zeros(num_state,1)].*impulse_response(ii,:),2),2);
        Cd(:,ii) = mod(sum([states ones(num_state,1)].*impulse_response(ii,:),2),2);
    end

    num_level = length(binary_data) / row;
    cost = inf(num_state, 1);
    cost(1) = 0;
    binary_data = reshape(binary_data,row,[])';
    idu = mod([0:num_state-1]'*2, num_state) + 1;
    idd = mod([0:num_state-1]'*2+1, num_state) + 1;
    previous = zeros(num_state, num_level);
    for ii = 1 : num_level
        costu = cost(idu) + sum(Cu ~= binary_data(ii,:), 2);
        costd = cost(idd) + sum(Cd ~= binary_data(ii,:), 2);
        [cost, idx] = min([costu costd], [], 2);
        idx = idx - 1;
        previous(:, ii) = idu.*(1 - idx) + idd.*(idx);
    end

    place = 1;
    decoded_data(num_level) = 0;
    for ii = num_level : -1 : 2
        place = previous(place, ii);
        decoded_data(ii-1) = (place > num_state/2);
    end
    decoded_data = decoded_data(1:(end-row));
end