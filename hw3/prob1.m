%% problem 1
binary_sequence1 = reshape(fliplr(de2bi([0:1]))',  1, []);
binary_sequence2 = reshape(fliplr(de2bi([0:3]))',  1, []);
binary_sequence3 = reshape(fliplr(de2bi([0:7]))',  1, []);
binary_sequence4 = reshape(fliplr(de2bi([0:15]))', 1, []);
d = 1;
symbol_sequence1  = symbol_mapper(binary_sequence1,  2, d, 'PAM', 'Binary');
symbol_sequence2  = symbol_mapper(binary_sequence2,  4, d, 'PAM', 'Binary');
symbol_sequence3  = symbol_mapper(binary_sequence3,  8, d, 'PAM', 'Binary');
symbol_sequence4  = symbol_mapper(binary_sequence4, 16, d, 'PAM', 'Binary');
symbol_sequence5  = symbol_mapper(binary_sequence1,  2, d, 'PSK', 'Binary');
symbol_sequence6  = symbol_mapper(binary_sequence2,  4, d, 'PSK', 'Binary');
symbol_sequence7  = symbol_mapper(binary_sequence3,  8, d, 'PSK', 'Binary');
symbol_sequence8  = symbol_mapper(binary_sequence4, 16, d, 'PSK', 'Binary');
symbol_sequence9  = symbol_mapper(binary_sequence2,  4, d, 'QAM', 'Binary');
symbol_sequence10 = symbol_mapper(binary_sequence4, 16, d, 'QAM', 'Binary');
symbol_sequence11 = symbol_mapper(binary_sequence1,  2, d, 'PAM', 'Gray');
symbol_sequence12 = symbol_mapper(binary_sequence2,  4, d, 'PAM', 'Gray');
symbol_sequence13 = symbol_mapper(binary_sequence3,  8, d, 'PAM', 'Gray');
symbol_sequence14 = symbol_mapper(binary_sequence4, 16, d, 'PAM', 'Gray');
symbol_sequence15 = symbol_mapper(binary_sequence1,  2, d, 'PSK', 'Gray');
symbol_sequence16 = symbol_mapper(binary_sequence2,  4, d, 'PSK', 'Gray');
symbol_sequence17 = symbol_mapper(binary_sequence3,  8, d, 'PSK', 'Gray');
symbol_sequence18 = symbol_mapper(binary_sequence4, 16, d, 'PSK', 'Gray');
symbol_sequence19 = symbol_mapper(binary_sequence2,  4, d, 'QAM', 'Gray');
symbol_sequence20 = symbol_mapper(binary_sequence4, 16, d, 'QAM', 'Gray');
%% testing error message
symbol_sequence21 = symbol_mapper(binary_sequence4, 14, d, 'QAM', 'Gray');

function symbol_sequence = symbol_mapper(binary_sequence, M, d, constellation, mapping)
    % – M: The number of points in the signal constellation.
        % M = 2,4,8,16 for 'PAM' & 'PSK'; M = 4,16 for 'QAM'.
    % – d: The minimum distance among the constellation.
    % – constellation: 'PAM', 'PSK' or 'QAM'.
    % – mapping: 'Binary' or 'Gray'.

    % For Gray code mapping of PAM, please assign 00 · · · 0 to the leftest constellation point.
    % Gray code mapping of PSK and QAM is shown in Figure 1 and Figure 2. If input M is
    % not the power of 2 or if input constellation is not defined, your function should be able
    % to throw error and display message.
    % Plot all the constellations, and show the bits of each symbol on the constellation points
    seq = [];
    mmax = 0;
    switch constellation
        case 'PAM'
            switch M
                case {2, 4, 8, 16}
                    mmax = 2;
                    seq = reshape(binary_sequence,log2(M),[])';
                    switch mapping
                        case 'Binary'
                            seq_len = length(seq);
                            C = cell(seq_len, 1);
                            for ii = 1:seq_len
                                C{ii,1} = (bi2de(fliplr(seq(ii,:)))*2-(M-1))*d;
                            end
                            symbol_sequence = [C{:}];
                        case 'Gray'
                            A2 = [1 0];
                            A4 = [2 1 3 0];
                            A8 = [5 4 2 3 6 7 1 0];
                            A16 = [0 1 3 2 7 6 4 5 15 14 12 13 8 9 11 10];
                            seq_len = length(seq);
                            C = cell(seq_len, 1);
                            for ii = 1:seq_len
                                idx = bi2de(fliplr(seq(ii,:)));
                                switch M
                                    case 2
                                        idx = A2(idx+1);
                                        C{ii,1} = idx*2-(M-1)*d;
                                    case 4
                                        idx = A4(idx+1);
                                        C{ii,1} = idx*2-(M-1)*d;
                                    case 8
                                        idx = A8(idx+1);
                                        C{ii,1} = idx*2-(M-1)*d;
                                    case 16
                                        idx = A16(idx+1);
                                        C{ii,1} = idx*2-(M-1)*d;
                                end
                            end
                            symbol_sequence = [C{:}];
                        otherwise
                            error("mapping must be 'Binary' or 'Gray'")
                    end
                otherwise
                    error("M must be 2, 4, 8, 16")
            end
        case 'PSK'
            switch M
                case {2, 4, 8, 16}
                    mmax = d;
                    seq = reshape(binary_sequence,log2(M),[])';
                    switch mapping
                        case 'Binary'
                            seq_len = length(seq);
                            C = cell(seq_len, 1);
                            for ii = 1:seq_len
                                idx = bi2de(fliplr(seq(ii,:)));
                                if M == 4
                                    C{ii,1} = cos(2*pi*idx/M+pi/4) + i*sin(2*pi*idx/M+pi/4);
                                else
                                    C{ii,1} = cos(2*pi*idx/M) + i*sin(2*pi*idx/M);
                                end
                            end
                            symbol_sequence = [C{:}];
                        case 'Gray'
                            A2 = [1 0];
                            A4 = [2 1 3 0];
                            A8 = [5 4 2 3 6 7 1 0];
                            A16 = [0 1 3 2 7 6 4 5 15 14 12 13 8 9 11 10];
                            seq_len = length(seq);
                            C = cell(seq_len, 1);
                            for ii = 1:seq_len
                                idx = bi2de(fliplr(seq(ii,:)));
                                switch M
                                    case 2
                                        idx = A2(idx+1);
                                        C{ii,1} = cos(2*pi*idx/M) + i*sin(2*pi*idx/M);
                                    case 4
                                        idx = A4(idx+1);
                                        C{ii,1} = cos(2*pi*idx/M+pi/4) + i*sin(2*pi*idx/M+pi/4);
                                    case 8
                                        idx = A8(idx+1);
                                        C{ii,1} = cos(2*pi*idx/M) + i*sin(2*pi*idx/M);
                                    case 16
                                        idx = A16(idx+1);
                                        C{ii,1} = cos(2*pi*idx/M) + i*sin(2*pi*idx/M);
                                end
                            end
                            symbol_sequence = [C{:}];
                        otherwise
                            error("mapping must be 'Binary' or 'Gray'")
                    end
                otherwise
                    error("M must be 2, 4, 8, 16")
            end
        case 'QAM'
            switch M
                case {4, 16}
                    mmax = (sqrt(M)-1)*d;
                    seq = reshape(binary_sequence,log2(M),[])';
                    switch mapping
                        case 'Binary'
                            seq_len = length(seq);
                            C = cell(seq_len, 1);
                            for ii = 1:seq_len
                                mid = log2(M)/2;
                                x_idx = bi2de(fliplr(seq(ii,1:mid)));
                                y_idx = bi2de(fliplr(seq(ii,mid+1:end)));
                                C{ii, 1} = ((x_idx*2-(sqrt(M)-1))*d) + i*((y_idx*2-(sqrt(M)-1))*d);
                            end
                            symbol_sequence = [C{:}];
                        case 'Gray'
                            A = [0 1 3 2]; 
                            seq_len = length(seq);
                            C = cell(seq_len, 1);
                            for ii = 1:seq_len
                                mid = log2(M)/2;
                                x_idx = bi2de(fliplr(seq(ii,1:mid)));
                                y_idx = bi2de(fliplr(seq(ii,mid+1:end)));
                                x_idx = A(x_idx+1);
                                y_idx = A(y_idx+1);
                                C{ii, 1} = ((x_idx*2-(sqrt(M)-1))*d) + i*((y_idx*2-(sqrt(M)-1))*d);
                            end
                            symbol_sequence = [C{:}];
                        otherwise
                            error("mapping must be 'Binary' or 'Gray'")
                    end
                otherwise
                    error("M must be 4, 16")
            end
        otherwise
            error("constellation must be 'PAM', 'PSK' or 'QAM'")
    end
    %% plot
    figure;
    seq_len = length(symbol_sequence);
    xx = roundn(real(symbol_sequence), -6);
    yy = roundn(imag(symbol_sequence), -6);
    plot(xx, yy, 'o', 'MarkerSize', 6, 'MarkerEdgeColor','b', 'MarkerFaceColor',[0.5,0.5,0.5]);
    xlabel('I channel');
    ylabel('Q channel');
    tit = strcat(num2str(M), '-', constellation, '-', mapping);
    A = cell(M,1);
    for ii = 1:M
        A{ii,1} = symbol_sequence(1,ii);
    end
    save(strcat('cell_prob1/', tit,'.mat'), 'A');
    title(tit);
    for ii = 1:seq_len
        text(xx(ii), yy(ii)+mmax/10, char(seq(ii,:)+48), 'Color', 'red', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
    saveas(gcf, strcat('image_prob1/', tit,'.png'));
end
