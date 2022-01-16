clear all;
close all;
oversampling_factor = 1000;
T_os = 1/oversampling_factor; % time spacing after oversampling
W = 50;
T = 1/(2*W);
pulse_duration = 1; % 1 sec
t_axis = (-pulse_duration/2 : T_os : pulse_duration/2 - T_os);
X1 = zeros(1, length(t_axis));
X1(1,(length(t_axis)/2-T/T_os):(length(t_axis)/2+T/T_os)) = 1;
X2 = sin(pi*t_axis/T)./(pi*t_axis/T); X2(1, find(isnan(X2))) = 1;
beta = 0.7;
X3 = X2.*cos(pi*beta*t_axis/T)./(1-4*beta^2*t_axis.^2/T^2); X3(1, find(isnan(X3))) = pi/4*sin(pi/2/beta)/(pi/2/beta);
% X4 = ( sin(pi*t_axis*(1-beta)/T) + 4*beta*t_axis/T.*cos(pi*t_axis/T*(1+beta)) ) ./ ( pi*t_axis/T.*(1-(4*beta*t_axis/T).^2) ) / sqrt(T);
% X4(1, find(isnan(X4))) = 1/sqrt(T)*(1+beta*(4/pi-1));
% X5 = X2 / sqrt(T);

[Wt, A1] = CTFT(X1, oversampling_factor, t_axis, 'Rect');
[Wt, A2] = CTFT(X2, oversampling_factor, t_axis, 'Sinc');
[Wt, A] = CTFT(X3, oversampling_factor, t_axis, 'Raised Cosine');

% [Wt, A] = CTFT(X4, oversampling_factor, t_axis, 'Root raised Cosine');
% [Wt, A] = CTFT(X5, oversampling_factor, t_axis, 'RSinc');
% X6 = conv(X5, X5)*T_os;
% [Wt, A] = CTFT(X6, oversampling_factor, (-pulse_duration+T_os:T_os:pulse_duration-T_os), 'SRSinc' );
% X7 = conv(X4, X4)*T_os;
% [Wt, A] = CTFT(X7, oversampling_factor, (-pulse_duration+T_os:T_os:pulse_duration-T_os), 'SRRaise' );

d = 2;
rng(19);
code = randi([0 1], 1, 20);
symbol_sequence = symbol_mapper_QPSK_gray(code, d);
Iphase = real(symbol_sequence);
Qphase = imag(symbol_sequence);

yI = pulse_shaper(Iphase, 'raised cosine');
t_axis2 = (0:T_os:length(yI)*T_os-T_os);
CTFT(yI, oversampling_factor, t_axis2, 'YI');

yQ = pulse_shaper(Qphase, 'raised cosine');
CTFT(yQ, oversampling_factor, t_axis2, 'YQ');

fc = 150;
t_axis2 = (0:T_os:length(yI)*T_os-T_os);
y = yI*sqrt(2).*cos(2*pi*fc*t_axis2) - yQ*sqrt(2).*sin(2*pi*fc*t_axis2);
CTFT(y, oversampling_factor, t_axis2, 'Y');

snrr = 25;
pd = makedist('Normal');
y = y + random(pd,1,length(y))*sqrt(sum(y.^2)/(10^(snrr/10))/length(y));
CTFT(y, oversampling_factor, t_axis2, 'Y + noise');
h = fir1(W, 0.9, 'low');

yIt = y*sqrt(2).*cos(2*pi*fc*t_axis2);
yIt = filter(h, 1, yIt);
CTFT(yIt, oversampling_factor, t_axis2, 'YI rcv');

yQt = -y*sqrt(2).*sin(2*pi*fc*t_axis2);
yQt = filter(h, 1, yQt);
CTFT(yQt, oversampling_factor, t_axis2, 'YQ rcv');

I_rcv = filter_sampling(yIt);
Q_rcv = filter_sampling(yQt);

function [W, A] = CTFT(x, fs, t_axis, tit)
    [col, row] = size(x);
    y = fft(x, row);

    N = row;
    w = 2*pi*(0:(N-1)) / N;
    w2 = fftshift(w);
    w3 = unwrap(w2 - 2*pi);

    y = y / fs;
    W = w3/pi*fs/2;
    A = abs(fftshift(y));
    figure;
    subplot(3,1,1);
    plot(t_axis, x);
    title(tit);
    
    subplot(3,1,2);
    plot(w3/pi*fs/2, abs(fftshift(y)));
    xlim([-500, 500]);
    title('CTFT Amplitude');
    xlabel('Freq (Hz)');
    ylabel('Amplitude');
    
    subplot(3,1,3);
    plot(w3/pi*fs/2, angle(fftshift(y)));
    xlim([-500, 500]);
    title('CTFT Phase');
    xlabel('Freq (Hz)');
    ylabel('Phase');
end

function y = pulse_shaper(x, pulse_shape, W)
    oversampling_factor = 1000;
    T_os = 1/oversampling_factor; % time spacing after oversampling
    W = 50;
    T = 1/(2*W);
    pulse_duration = 1; % 1 sec
    t_axis = (-pulse_duration/2 : T_os : pulse_duration/2 - T_os);
    Rect = zeros(1, length(t_axis));
    Rect(1,(length(t_axis)/2-T/T_os):(length(t_axis)/2+T/T_os)) = 1;
    RSinc = sin(pi*t_axis/T)./(pi*t_axis/T) / sqrt(T);
    RSinc(1, find(isnan(RSinc))) = 1 / sqrt(T);
    beta = 0;
    RRaised = ( sin(pi*t_axis*(1-beta)/T) + 4*beta*t_axis/T.*cos(pi*t_axis/T*(1+beta)) ) ./ ( pi*t_axis/T.*(1-(4*beta*t_axis/T).^2) ) / sqrt(T);
    RRaised(1, find(isnan(RRaised))) = 1/sqrt(T)*(1+beta*(4/pi-1));

    switch pulse_shape
    case 'sinc'
        Pul = RSinc;
    case 'raised cosine'
        Pul = RRaised;
    otherwise
        error("Only 'sinc' and 'raised cosine' are valid for pulse shape.")
    end

    save('Pulse.mat', 'Pul');
    O = zeros(1, (T/T_os)*(length(x)-1)+length(Pul));
    O(1:length(Pul)) = Pul;
    y = zeros(1, length(O));
    for ii = 1 : length(x)
        y = y + x(1,ii)*circshift(O,(T/T_os)*(ii-1));
    end

end

function symbol_sequence = symbol_mapper_QPSK_gray(binary_sequence, d)
    seq = [];
    M = 4;
    seq = reshape(binary_sequence,log2(M),[])';
    A4 = [2 1 3 0];
    seq_len = length(seq);
    C = cell(seq_len, 1);
    for ii = 1:seq_len
        idx = bi2de(fliplr(seq(ii,:)));
        idx = A4(idx+1);
        C{ii,1} = cos(2*pi*idx/M+pi/4) + i*sin(2*pi*idx/M+pi/4);
    end
    symbol_sequence = [C{:}];
end

function y = filter_sampling(in)
    oversampling_factor = 100000;
    T_os = 1/oversampling_factor; % time spacing after oversampling
    C = {};
    A = load('Pulse.mat');
    Pul = A.Pul;
    RE = conv(in, Pul)*T_os;
    ii = 1;
    while ii*length(Pul) <= length(RE)
        C{ii} = RE(1,ii*length(Pul));
        ii = ii + 1;
    end
    y = [C{:}];
end