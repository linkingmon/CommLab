% prolbem 1 

R = [10 100 1000];                                  % different simulation times
N = 1000;                                           % number of bits
sym = load('../hw3/cell_prob1/16-QAM-Gray.mat');    % load constellation
sym = [sym.A{:}];
EboN0_list = 10.^([0:30]/10);
ber_res = zeros(length(R), length(EboN0_list));
ser_res = zeros(length(R), length(EboN0_list));

rng(2);
for ii = 1 : length(R)
    for jj = 1 : length(EboN0_list)
        EboN0 = EboN0_list(jj);
        message = randi([0 1], R(ii), N);
        symbol_map = reshape(bi2de(fliplr(reshape(message', 4, [])'))' + 1, N/4, [])';
        symbol = sym(symbol_map);
        Es = sum(abs(symbol).^2, 2) / (N/4);
        noise_var = Es / (4*EboN0);

        pd = makedist('Normal');
        noise = random(pd, R(ii), N / 4).*sqrt(noise_var/2) + i * random(pd, R(ii) ,N / 4).*sqrt(noise_var/2);
        symbol_rec = symbol + noise;
        [symbol_demap, message_rec] = demap(symbol_rec);
        
        ber_res(ii, jj) = sum(message ~= message_rec, 'all') / numel(message);
        ser_res(ii, jj) = sum(symbol_map ~= symbol_demap, 'all') / numel(message) * 4;

        % figure;
        % xx = roundn(real(symbol_rec), -6);
        % yy = roundn(imag(symbol_rec), -6);
        % plot(xx, yy, 'o', 'MarkerSize', 6, 'MarkerEdgeColor','b', 'MarkerFaceColor',[0.5,0.5,0.5]);
        % hold on;
        % xx = roundn(real(symbol), -6);
        % yy = roundn(imag(symbol), -6);
        % plot(xx, yy, 'o', 'MarkerSize', 8, 'MarkerEdgeColor','r', 'MarkerFaceColor',[1,0,0]);
    end
end

save('problem1_result.mat', 'ber_res', 'ser_res');

% 16-QAM demapper
function [sym_out, bin_out] = demap(sym_in)
    persistent sym
        if isempty(sym)
            sym = load('../hw3/cell_prob1/16-QAM-Gray.mat');
            sym = [sym.A{:}];
        end
    [row, sym_len] = size(sym_in);
    C = zeros(row, sym_len);
    for ii = 1 : sym_len
        [minvalue, index_of_min] = min(abs(sym_in(:, ii) - sym), [], 2);
        C(:,ii) = index_of_min - 1;
    end
    sym_out = C + 1;
    bin_out = reshape(fliplr(de2bi(C',4))',sym_len*4,[])';
end