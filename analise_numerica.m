%% Fast Fourier Transform - Análise Numérica
% Instituto Federal Fluminense - Engenharia de Controle e Automação
% Processamento de Sinais (2022.1)
% Prof.: Alexandre C. Leite
% Aluno: Kaique Guimarães Cerqueira
clear, close all

%% Dado o sinal de entrada u(t):
% Frequência de amostragem (em Hz) - ao menos 2 vezes maior que a maior
%    30 Hz                            componente em frequência do sinal
Fs = 30;
Ts = 1/Fs;
% Quantidade de pontos:
L = 2^12;
t = (0:(L-1))*Ts; % Vetor de tempo
%% Geração de dados: 
% % Foi utilizado o seguinte código inicialmente para gerar os dados. Para
% % gerar novas amostras ou dados com mais pontos, descomentar esta parte.
% % Sinal de entrada: 
% % 10Hz(sin defasado em 45º) + 3Hz(cos c maior pot.) + 14Hz(sin menor pot.)
% u_t_clean = sin(2*pi*10*t + pi/2) + 1.5*cos(2*pi*3*t) + 0.75*sin(2*pi*14*t);
% % Adicionando ruído de distribuição gaussiana:
% u_t = u_t_clean + 2.5*randn(size(t));

% Sinal carregado (4096 pontos):
load('signal_data.mat');
figure()
plot(t, u_t, 'r', 'LineWidth', 1.75);
hold on
plot(t, u_t_clean, '-k', 'LineWidth', 1);
xlabel("Tempo (s)")
ylabel("Amplitude")
legend('Sinal com ruído', 'Sinal Limpo')
axis([0 2*pi 1.5*min(u_t) 1.5*max(u_t)]);

%% FFTs
% do MATLAB
tic;
Y_matlab = fft(u_t);
toc;
Y_matlab = fftshift(Y_matlab);  % Centralizando a frequência em 0

% FFT iterativa em C
Y_c = myfft(u_t);
Y_c = fftshift(Y_c);            % Centralizando a frequência em 0

% Power Spectrum function (módulo normalizado)
PSD = @(x) abs(x).^2/length(x);

powershift_matlab = PSD(Y_matlab);
powershift_c = PSD(Y_c);        

fshift = (-L/2:L/2-1)*(Fs/L); % zero-centered frequency range

% Plots das FFTs
figure()
plot(fshift,powershift_matlab, 'r', 'LineWidth', 2.5)
% title('FFT do MATLAB')
xlabel('Frequência (Hz)')
hold on
% figure()
plot(fshift,powershift_c, '--k', 'LineWidth', 1)
% title('FFT implementada em C')
legend('FFT do MATLAB', 'FFT implementada em C', 'Location','best')
xlabel('Frequência (Hz)')

%% Análise Numérica
% Métricas: Erro médio quadrático
emq = immse(Y_matlab,Y_c);
sprintf('Erro medio quadrático: %f', emq)

% Erro máximo
err_max = max(abs(Y_matlab - Y_c));
sprintf('Erro máximo: %f', err_max)

% Plot do erro em módulo
figure()
plot(fshift, abs(Y_matlab - Y_c))
ylabel("$|Y_i - \hat{Y}_i|$", 'Interpreter','latex')
xlabel("Frequência (Hz)")
