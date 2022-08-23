function output = myfft(samples)
%% The Fast Fourier Transform Algorithm
% Kaique Guimarães Cerqueira
% Objectives: 
%   [x] To read and write values in double format into a .bin file
%   [x] To implement iteractive FFT algorithm
%   [ ] To implement Cooley-Tukey FFT algorithm in C (recursively)
%   [ ] To be able to compare results both in C and MATLAB (maybe other libs)
%       - Parameters:
%           Time
%           Precision

%% Writing
fileID = fopen('sample.bin', 'w');
[~,~,~,encoding] = fopen(fileID); % Checking encoding system
fwrite(fileID, samples,"double",'n'); % utilizando o formato nativo da máquina
fclose(fileID);
%% Call program in C
% Iterative Program
system('itfft.exe');
% Recursive Program (to be implemented)

%% Reading
fileID = fopen('output.bin', 'r');
[~,~,~,encoding] = fopen(fileID); % Checking encoding system
samples = fread(fileID, 'double'); % Reading binary and using as double
fclose(fileID);

tamanho = (numel(samples)/2); % Definindo o tamanho do vetor resposta
output = samples(1:tamanho)'*0; % Inicializando vetor
for k = 1:tamanho
    output(k) = samples(k*2-1) + 1i*samples(k*2);
end

end