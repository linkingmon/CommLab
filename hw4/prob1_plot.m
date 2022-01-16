% plot the BER-(Eb/N0) curve

% load ber data
T = load('problem1_result.mat');
EboN0_list = [0:30];        % dB
ber_res = T.ber_res;
ser_res = T.ser_res;
SNR = 10.^(EboN0_list/10)*4;
ser_formula = qfunc(sqrt(SNR/5));

% plot BER - (Eb/N0)
figure;
plot(EboN0_list, ber_res(1,:), '-', 'color', 'g', 'LineWidth', 1);
hold on;
plot(EboN0_list, ber_res(2,:), '-', 'color', 'r', 'LineWidth', 1);
plot(EboN0_list, ber_res(3,:), '-', 'color', 'b', 'LineWidth', 1);
legend('10 times', '100 times', '1000 times');
set(gca, 'YScale', 'log')
title('BER - (Eb/N0)');
ylabel('BER');
xlabel('(Eb/N0) [dB]');
hold off;

% plot SER - (Eb/N0)
figure;
plot(EboN0_list, ser_res(1,:), '-', 'color', 'g', 'LineWidth', 1);
hold on;
plot(EboN0_list, ser_res(2,:), '-', 'color', 'r', 'LineWidth', 1);
plot(EboN0_list, ser_res(3,:), '-', 'color', 'b', 'LineWidth', 1);
legend('10 times', '100 times', '1000 times');
set(gca, 'YScale', 'log')
title('SER - (Eb/N0)');
ylabel('SER');
xlabel('(Eb/N0) [dB]');
hold off;

% plot SER - (Eb/N0)
figure;
plot(EboN0_list, ser_res(3,:), '-', 'color', 'b', 'LineWidth', 1);
hold on;
plot(EboN0_list, 3*ser_formula-9/4*ser_formula.^2, '-', 'color', 'r', 'LineWidth', 1);
plot(EboN0_list, ser_formula, '-', 'color', 'g', 'LineWidth', 1);
plot(EboN0_list, 15*ser_formula, '-', 'color', 'y', 'LineWidth', 1);
legend('simulated', 'theoretical', 'lower bound', 'upper bound');
set(gca, 'YScale', 'log')
title('SER - (Eb/N0)');
xlim([0,15])
ylabel('SER');
xlabel('(Eb/N0) [dB]');
hold off;
