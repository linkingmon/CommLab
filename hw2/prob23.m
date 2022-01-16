%% problem 2

% assume amp is in [-1, 1]
[x, fs] = audioread('handel.ogg');
xmax = 1;
bit = 4;
level = 2^bit;
xt = quantizer_L_level(x, xmax, level)';

% (e) PCM
y = pcm_enc(xt, 4);
length_y = length(y)

% (f) huffman dict
delta = 2 * xmax / level;
symbols = [-(level-1)*delta/2:delta:(level-1)*delta/2];
p = histc(xt, symbols);
p = p / sum(p);
dict = huffman_dict(symbols, p);
dict2 = huffmandict(symbols, p);

% (g) huffman encode
y = huffmanenco(xt, dict);
length_y = length(y);
y1 = huffman_enc(xt, dict);
length_y = length(y1);
check_enc = sum(y ~= y1)            % 0 for correct

% (h) huffman decode
xd = huffmandeco(y, dict);
length_x = length(xd);
xd2 = huffman_dec(y, dict);
length_x = length(xd2);
check_dec = sum(xd ~= xd2)          % 0 for correct

% clear all;

function y = quantizer_L_level(x, xmax, level)
    delta = 2 * xmax / level;
    partition = [-xmax:delta:xmax];
    codebook = [0,-(level-1)*delta/2:delta:(level-1)*delta/2,0];
    [I, y] = quantiz(x,partition,codebook); 
end

function y = pcm_enc(x, numBits)
    level = 2^numBits;
    delta = 2 / level;
    x = (x + (level-1)*delta/2) / delta + 1;
    y = reshape(de2bi(x)',1,[]);
end

function dict = huffman_dict(symbols, p)
    dict = cell(length(symbols), 2);
    for i = 1 : length(symbols)
        dict{i,1} = symbols(i);
    end
    % node is a structure array
    for i = 1 : length(symbols)
        s.array = [i];
        s.prob = p(i);
        node(i) = s;
    end
    nodeb = node
    while length(nodeb) > 1
        min = 1.1;
        minmin = 1.1;
        minp = -1;
        minminp = -1;
        for i = 1:length(nodeb)
            if nodeb(i).prob < minmin
                min = minmin;
                minp = minminp;
                minmin = nodeb(i).prob;
                minminp = i;
                continue
            end
            if nodeb(i).prob < min
                min = nodeb(i).prob;
                minp = i;
            end
        end
        if minp > minminp
            temp = minp;
            minp = minminp;
            minminp = temp;
        end
        for i = 1:length(nodeb(minp).array)
            dict{nodeb(minp).array(i),2} = [0 dict{nodeb(minp).array(i),2}];
        end
        for i = 1:length(nodeb(minminp).array)
            dict{nodeb(minminp).array(i),2} = [1 dict{nodeb(minminp).array(i),2}];
        end
        s.array = [nodeb(minp).array nodeb(minminp).array];
        s.prob = nodeb(minp).prob + nodeb(minminp).prob;
        nodeb = [nodeb(1:minp-1), nodeb(minp+1:minminp-1), nodeb(minminp+1:end)];
        nodeb = [nodeb s];
        nodeb_length = length(nodeb);
    end
end

function y = huffman_enc(x, dict)
    C = cell(length(x),1);
    for i = 1 : length(x)
        C{i, 1} = dict{find( [dict{:,1}] == x(i)),2};
    end
    y = [C{:,1}]';
end

function index = code2idx(codeword)
    index = 1;
    for i = 1: length(codeword)
        index = index*2;
        if codeword(i) == 0
            index = index +1;
        end
    end
end

function tree = constructHuffmanTree(dict)
    for i = 1:length(dict)
        leaf = code2idx(dict{i,2});
        tree{leaf} = dict{i,1};
    end
end

function x_dec = huffman_dec(y, dict)
    x_dec{393586} = [];
    tree = constructHuffmanTree(dict);
    codepos = 1; % position of the codeword
    codelen = 1; % length of the codeword
    x_idx = 1;
    while codepos < length(y)
        while 1
            idx = code2idx(y(codepos:codepos+codelen-1)); % the position of the tree
            if isempty(tree{idx})
                codelen = codelen +1;
            else
                x_dec{x_idx} = tree{idx};
                x_idx = x_idx+1;
                codepos = codepos + codelen;
                codelen =1;
                break;
            end
        end
    end
end

% function x_dec = huffman_dec(y, dict)
%     s.array = [1:length(dict)];
%     s.level = 1;
%     % 0 : left 1 : right
%     s = split_tree(s, dict);
%     st = s;
%     D = cell(length(y),1);
%     cnt = 1;
%     for i = 1 : length(y)
%         i
%         if y(i) == 1
%             st = st.right;
%             if length(st.array) == 1
%                 D{cnt, 1} = st.array(1);
%                 st = s;
%                 cnt = cnt + 1;
%             end
%         else
%             st = st.left;
%             if length(st.array) == 1
%                 D{cnt, 1} = st.array(1);
%                 st = s;
%                 cnt = cnt + 1;
%             end
%         end
%     end
%     x_dec = [D{:, 1}];
%     for i = 1 : length(x_dec)
%         x_dec(i) = dict{x_dec(i),1};
%     end
% end

% function st = split_tree(s, dict)
%     if length(s.array) <= 1
%         st = s;
%         return
%     end
%     x.array = [];
%     y.array = [];
%     for i = 1 : length(s.array)
%         if dict{s.array(i),2}(s.level) == 0
%             x.array = [x.array s.array(i)];
%         else
%             y.array = [y.array s.array(i)];
%         end
%     end
%     x.level = s.level + 1;
%     y.level = s.level + 1;
%     st.array = s.array;
%     st.left = split_tree(x, dict);
%     st.right = split_tree(y, dict);
% end

% function x_dec = huffman_dec(y, dict)
%     x_dec = []
%     yt = y;
%     pos = 1;
%     dict2.symbol = dict(:,1);
%     dict2.code = dict(:,2);
%     i = 1
%     while( ~isempty(yt) )
%         temp = yt(pos);                                 % Get first bit(char). 
%         dictb = dict2;                                   % Get a backup of the input dictionary.
%         while(1)
%             dictb = found_match( temp, pos, dictb);     % Get a sub dictionary.
%             if ( length(dictb.code) ~= 1 )              % Until one codeword left.
%                 pos = pos + 1;                          % Match second bit(char).
%                 temp = yt(pos);                         % Get first char 
%             else
%                 pos = 1;                                % Reset position.
%                 yt = yt(length(dictb.code{1})+1:end);   % Update the input signal.
%                 break;
%             end        
%         end
%         i = i + 1
%         % x_dec = [x_dec dictb.symbol]                     % Append char to decoded signal.
%         % size(x_dec)
%     end
% end

% function dictb = found_match( code, pos, dict )
%     dictb.symbol={}; dictb.code={};                     % Create a dictionary structure.
%     j = 1;                                              % Iterator.
%     for i = 1:length(dict.code)                         % For each code in dictionary.
%         if ( dict.code{i}(pos) == code )                % If inpute code matches 
%              dictb.symbol(j) = dict.symbol(i);          % Get the symbol that matches.
%              dictb.code(j) = dict.code(i);              % Get the code of the matched symbol.
%              j = j + 1;                                 % Prepare for a next symbol.
%         end
%     end
% end