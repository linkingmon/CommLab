%% problem 4
close all;
clear all;
binary_sequence = reshape(fliplr(de2bi(randi(16, 1, 2400)-1,4))',  1, []);
d = 1;
berr_all = {};
serr_all = {};
for SNR=0:1:20   
    cnt = 0;
    tit = cell(1,50);
    berr = cell(1,50);
    serr = cell(1,50);
    % [tit{1}, berr{1}, serr{1}] = Run(binary_sequence,  2, d, 'PAM', 'Gray', 'MD', SNR);
    % [tit{2}, berr{2}, serr{2}] = Run(binary_sequence,  4, d, 'PAM', 'Gray', 'MD', SNR);
    % [tit{3}, berr{3}, serr{3}] = Run(binary_sequence,  8, d, 'PAM', 'Gray', 'MD', SNR);
    % [tit{4}, berr{4}, serr{4}] = Run(binary_sequence, 16, d, 'PAM', 'Gray', 'MD', SNR);
    % [tit{5}, berr{5}, serr{5}] = Run(binary_sequence,  2, d, 'PSK', 'Gray', 'MD', SNR);
    % [tit{6}, berr{6}, serr{6}] = Run(binary_sequence,  4, d, 'PSK', 'Gray', 'MD', SNR);
    % [tit{7}, berr{7}, serr{7}] = Run(binary_sequence,  8, d, 'PSK', 'Gray', 'MD', SNR);
    % [tit{8}, berr{8}, serr{8}] = Run(binary_sequence, 16, d, 'PSK', 'Gray', 'MD', SNR);
    [tit{9}, berr{9}, serr{9}] = Run(binary_sequence,  4, d, 'QAM', 'Gray', 'MD', SNR);
    [tit{10}, berr{10}, serr{10}] = Run(binary_sequence, 16, d, 'QAM', 'Gray', 'MD', SNR);
    % [tit{11}, berr{11}, serr{11}] = Run(binary_sequence,  2, d, 'PAM', 'Binary', 'MD', SNR);
    % [tit{12}, berr{12}, serr{12}] = Run(binary_sequence,  4, d, 'PAM', 'Binary', 'MD', SNR);
    % [tit{13}, berr{13}, serr{13}] = Run(binary_sequence,  8, d, 'PAM', 'Binary', 'MD', SNR);
    % [tit{14}, berr{14}, serr{14}] = Run(binary_sequence, 16, d, 'PAM', 'Binary', 'MD', SNR);
    % [tit{15}, berr{15}, serr{15}] = Run(binary_sequence,  2, d, 'PSK', 'Binary', 'MD', SNR);
    % [tit{16}, berr{16}, serr{16}] = Run(binary_sequence,  4, d, 'PSK', 'Binary', 'MD', SNR);
    % [tit{17}, berr{17}, serr{17}] = Run(binary_sequence,  8, d, 'PSK', 'Binary', 'MD', SNR);
    % [tit{18}, berr{18}, serr{18}] = Run(binary_sequence, 16, d, 'PSK', 'Binary', 'MD', SNR);
    % [tit{19}, berr{19}, serr{19}] = Run(binary_sequence,  4, d, 'QAM', 'Binary', 'MD', SNR);
    % [tit{20}, berr{20}, serr{20}] = Run(binary_sequence, 16, d, 'QAM', 'Binary', 'MD', SNR);
    % [tit{21}, berr{21}, serr{21}] = Run(binary_sequence, 16, d, 'PAM', 'Gray', 'MAP', SNR);
    % [tit{22}, berr{22}, serr{22}] = Run(binary_sequence, 16, d, 'PSK', 'Gray', 'MAP', SNR);
    % [tit{23}, berr{23}, serr{23}] = Run(binary_sequence, 16, d, 'QAM', 'Gray', 'MAP', SNR);
    % [tit{24}, berr{24}, serr{24}] = Run(binary_sequence, 16, d, 'PAM', 'Gray', 'ML', SNR);
    % [tit{25}, berr{25}, serr{25}] = Run(binary_sequence, 16, d, 'PSK', 'Gray', 'ML', SNR);
    % [tit{26}, berr{26}, serr{26}] = Run(binary_sequence, 16, d, 'QAM', 'Gray', 'ML', SNR);
    
    % figure;
    berr = [berr{:}];
    serr = [serr{:}];
    berr_all{SNR+1} = berr;
    serr_all{SNR+1} = serr;
    tit = tit(~cellfun('isempty',tit));
    % C = categorical(tit);
    % h = histogram('Categories', C, 'BinCounts', berr, 'FaceColor', 'b', 'FaceAlpha', 1);
    % title(strcat('SNR = ',{' '},num2str(SNR),'dB'));
    % ylabel('Bit error rate');
    % xlabel('constellation & mapping');
end
n = length(tit);
A = reshape([berr_all{1,:}],n,[]);
figure;
plot([0:1:20], A(1,:));
xlabel('SNR');
ylabel('bit error rate');
title('BER - SNR')
set(gca, 'YScale', 'log')
hold on;
for ii = 2:n
    plot([0:1:20], A(ii,:));
end
legend(tit,'Location','northeast')
hold off;
saveas(gcf, strcat('image_prob4/SNRplot/', [tit{:}], '_BER','.png'));

A = reshape([serr_all{1,:}],n,[]);
figure;
plot([0:1:20], A(1,:));
xlabel('SNR');
ylabel('symbol error rate');
title('SER - SNR')
set(gca, 'YScale', 'log')
hold on;
for ii = 2:n
    plot([0:1:20], A(ii,:));
end
legend(tit,'Location','northeast')
hold off;
saveas(gcf, strcat('image_prob4/SNRplot/', [tit{:}],'_SER','.png'));


function [tit, berr, serr] = Run(binary_sequence, M, d, constellation, mapping, decision_rule, SNR)
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
    binary_sequence_in_symbol = bi2de(reshape(binary_sequence, log2(M), [])');
    binary_sequence_rec_in_symbol = bi2de(reshape(binary_sequence_rec, log2(M), [])');
    serr = sum(binary_sequence_in_symbol ~= binary_sequence_rec_in_symbol) / ( binary_seq_len/log2(M) );
    fprintf("%16s Bit error rate   %f\n", tit, berr);
    fprintf("%16s Symbol error rate   %f\n", tit, serr);
    tit = strcat(tit, '-', decision_rule, '-', int2str(SNR));

    % figure;
    % seq_len = length(symbol_sequence);
    % xx = roundn(real(symbol_sequence_rec), -6);
    % yy = roundn(imag(symbol_sequence_rec), -6);
    % plot(xx, yy, 'o', 'MarkerSize', 6, 'MarkerEdgeColor','b', 'MarkerFaceColor',[0.5,0.5,0.5]);
    % hold on;
    % xx = roundn(real(symbol_sequence), -6);
    % yy = roundn(imag(symbol_sequence), -6);
    % plot(xx, yy, 'o', 'MarkerSize', 8, 'MarkerEdgeColor','r', 'MarkerFaceColor',[1,0,0]);
    % xlabel('I channel');
    % ylabel('Q channel');
    % A = cell(M,1);
    % for ii = 1:M
    %     A{ii,1} = symbol_sequence(1,ii);
    % end
    % title(tit);
    % saveas(gcf, strcat('image_prob4/', tit,'.png'));
    % hold off;
end

function binary_sequence = symbol_demapper(symbol_sequence, M, d, constellation, mapping, decision_rule)
    % where decision_rule is 'MD'.
    % (Bonus) decision_rule: 'MAP' and 'ML'.
    
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
                % abs(B - symbol_sequence(ii))
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
                % abs(B - symbol_sequence(ii))
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