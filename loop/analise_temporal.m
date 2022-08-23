%% Fast Fourier Transform - Análise Temporal
% Instituto Federal Fluminense - Engenharia de Controle e Automação
% Processamento de Sinais (2022.1)
% Prof.: Alexandre C. Leite
% Aluno: Kaique Guimarães Cerqueira
clear, close all

% Computational samples
profiling = 1:10000;
time_itfft = profiling*0;
% Run once in 19/08/2022. Exp time 9h
% for i = profiling
% %     Run iteractive FFT 1000 times (in C code, it is already being run
% %     100 times)
%     [~,cmdoutput] = system("loop_itfft.exe");
%     time_itfft(i) = str2num(cmdoutput);
%     i
% end
% save('itfft_profiling', 'time_itfft');

% Loading data
load("itfft_profiling.mat");
%% Data analysis
% Calculating mean value for each batch
time_itfft = time_itfft/100;
% Mean
t_media = mean(time_itfft)
% Standard deviation
t_desviopadrao = std(time_itfft)
% Variance
t_var = var(time_itfft)
% Minimum and maximum ocurrencies
t_min = min(time_itfft)
t_max = max(time_itfft)

%% Plots: histogram
[counts, binCenters] = hist(time_itfft, 1000); %#ok<HIST> 
% Plot da probabilidade relativa temporal de cada FFT 
figure()
bar(binCenters, counts/length(counts))
ylabel('Probabilidade')
xlabel('Tempo (s)')
pause;
% Ajuste de eixo
axis([0.02 0.022 0 1])

%% Gaussian Analysis - made using cftool
[xData, yData] = prepareCurveData( binCenters, counts );
% 
% % Set up fittype and options.
% ft = fittype( 'gauss8' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
% opts.Normalize = 'on';
% opts.StartPoint = [771 -0.115989423886605 0.0180710324125379 314.482049826262 -0.102139940437458 0.0213524690317421 110.340674373789 -0.140226019922612 0.0271184646489327 22.6236970080996 -0.0467420066408708 0.0419209782294982 20 -1.69483053708934 0.0687670604311507 18.054104479099 0.00865592715571637 0.0400866352340347 8.82929809552389 0.0952151987128842 0.058262536079014 6.99986319346555 0.28910796700094 0.0881070458364841];
% 
% % Generate fit model (Gaussian w/ 8 terms)
% [gaussian_8, gof] = fit( xData, yData, ft, opts );
% Gaussian_8 coeffs
       a1 =       748.8;  
       b1 =     0.02086;  
       c1 =   0.0001824;  
       a2 =       32.27;  
       b2 =     0.02087;  
       c2 =   0.0004331;  
       a3 =       74.43;  
       b3 =     0.02055;  
       c3 =   0.0001032;  
       a4 =       9.984;  
       b4 =     0.02162;  
       c4 =   0.0007228;  
       a5 =       8.714;  
       b5 =    0.006366;  
       c5 =   0.0006109;  
       a6 =        8.03;  
       b6 =     0.02197;  
       c6 =   0.0001355;  
       a7 =       2.578;  
       b7 =      0.0234;  
       c7 =    0.001759;  
       a8 =       2.073;  
       b8 =     0.02444; 
       c8 =    0.000431;

gaussian_8 = @(x) a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2) + a3*exp(-((x-b3)/c3).^2) + a4*exp(-((x-b4)/c4).^2) + a5*exp(-((x-b5)/c5).^2) + a6*exp(-((x-b6)/c6).^2) + a7*exp(-((x-b7)/c7).^2) + a8*exp(-((x-b8)/c8).^2)
% Calculate integral from fitted curve
int = integral(gaussian_8, 0, inf);
% Normalized curve
gaussian_8_norm = @(x) gaussian_8(x)/int;
% Now we can calculate the probability of the events inside the boundaries
% delimited by the integrals by defining its limits:
fastest = integral(gaussian_8_norm, 0, 0.02); % Prob of being less than 0.02s
medio = integral(gaussian_8_norm, 0.02, 0.022); % Prob of being in middle
slowest = integral(gaussian_8_norm, 0.025, inf); % Prob of being more than 0.025s
sprintf(['Probability of being less than 0.02s: %f %%\nProbability of ' ...
    'being more than 0.025s: %f %%'], fastest*100, slowest*100)

% Plot fit with data.
figure( 'Name', 'gaussian_8' );
plot(xData, gaussian_8(xData), 'LineWidth', 2.5);
hold on
plot(xData, yData, '.r', 'MarkerSize', 8);
legend('Gaussian curve model w/ 8 terms', 'Frequency of observations', 'Location', 'northeast');
title('FFT Profiling: Duração do algoritmo')
% Label axes
xlabel( 'tempo (s)', 'Interpreter', 'none' );
ylabel( 'Nº de eventos', 'Interpreter', 'none' );
axis([xData(1) xData(end) 0 800])
grid on

clear opts ft xData yData
%% Analise temporal: FFT do Matlab
time_matlab = profiling*0;
% Run a million times
% for i = profiling
%     % Loading benchmark signal
%     load('..\signal_data.mat')
%     tic
%     for j = 1:100 % Run FFT matlab function 100 times
%         fft(u_t);
%     end
%     time_matlab(i) = toc;
%     clear u_t u_t_clean;
% end
% save('Mfft_profiling', 'time_matlab');

load('Mfft_profiling', 'time_matlab');
time_matlab = time_matlab / 100;
mat_mean = mean(time_matlab);
%% Analisar DFT