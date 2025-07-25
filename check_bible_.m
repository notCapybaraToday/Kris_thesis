clc; clear; close all;

% Input values
L = 0.7;
R1 = 0.0021;
R2 = 0.0235;
t = 25;  % Celsius
mu = 1.825e-5;
nfreq = 3000;
f = linspace(20, 1000, nfreq);

% Call the function for book implementation
[Z_in_lossy, Z_in_lossless, f] = compute_input_impedance(f, L, R1, R2, t, mu);
% Call the function for article implementation
[Z_in_lossy1, Z_in_lossless1, f1] = input_impedance_article(f, L, R1, R2, t, mu);
% Call the function to read CSV from openwind
[freq, Z] = open_wind_csv('Impedance_25degC (1)_05Hz.csv');




% ----- comparison with losses 
figure;
subplot(2,1,1);
plot(f, db(Z_in_lossy1), 'LineWidth', 2.5); hold on;
plot(f, db(Z_in_lossy), 'LineWidth', 1.5);
plot(freq,db(Z),'--', 'LineWidth', 3);
legend('With losses -a', 'With losses -b', 'With losses - open wind');
xlabel('Frequency (Hz)');
ylabel('|Z_{in}/Z_{c}| (dB)');
title('Normalized Input Impedance Magnitude Comparison (with losses)');
grid on;

subplot(2,1,2);
plot(f, angle(Z_in_lossy1), 'LineWidth', 2.5); hold on;
plot(f, angle(Z_in_lossy), 'LineWidth', 1.5);
plot(freq,angle(Z), '--','LineWidth', 3);
legend('With losses -a', 'With losses -b', 'With losses - open wind');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
title('Input Impedance Phase Comparison (with losses)');
grid on;

% ------- comparison NO losses
figure;
subplot(2,1,1);
plot(f, db(Z_in_lossless1), 'LineWidth', 2.5); hold on;
plot(f, db(Z_in_lossless), 'LineWidth', 1.5);
plot(freq,db(Z),'--', 'LineWidth', 3);
legend('Lossless -a', 'Lossless -b', 'With losses - open wind');
xlabel('Frequency (Hz)');
ylabel('|Z_{in}| (dB)');
title('Normalized Input Impedance Magnitude Comparison (lossless)');
grid on;

subplot(2,1,2);
plot(f, angle(Z_in_lossless1),'LineWidth', 2.5); hold on;
plot(f, angle(Z_in_lossless), 'LineWidth', 1.5);
plot(freq,angle(Z),'--', 'LineWidth', 3);
legend('Lossless -a', 'Lossless -b', 'With losses - open wind');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
title('Input Impedance Phase Comparison (lossless)');
grid on;


%------ error size
compare_resonances(f, Z_in_lossy1, f, Z_in_lossless, freq, Z, 4);

