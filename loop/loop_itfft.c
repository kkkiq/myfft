/*
        Engenharia de Controle e Automação - Instituto Federal Fluminense
        Disciplina: Processamento de Sinais (2022.1)
        Prof.: Alexandre C. Leite
        Aluno: Kaique Guimarães Cerqueira

        @summary
        Código para rodar em batches a implementação iterativa da FFT, a
        fim de computar métricas temporais.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

int main() {
    // ======================================================================
    //                      Definindo tamanho do arquivo
    // ======================================================================
    FILE* pSize = fopen("benchmark.bin", "rb");
    fseek(pSize, 0L, SEEK_END); // Procura o final do arquivo
    long file_size = ftell(pSize); // Tamanho (em bytes)
    // printf("O arquivo possui %ld bytes.\n", file_size);
    fclose(pSize);
    // Número de amostras:
    const long raw_sample_size = (long)(file_size/(sizeof(double)));
    long sample_size = 1;
    // Definição do vetor com tamanho de 2^N amostras
    while (sample_size < raw_sample_size)
    {
        sample_size *= 2;
    }
    
    // ======================================================================
    //                            Start Loop
    // ======================================================================
    // Definição das variáveis de tempo:
    clock_t start = 0, end = 0;
    double timer = 0;

    for (int j = 0; j < 1000; j++)
    {
        // ======================================================================
        //                         Leitura do arquivo
        // ======================================================================
        char buffer[8]; // Buffer das amostras (64 bits cada)
        double converted_data; // Dado após ser convertido em double
        complex samples[sample_size]; // Array de amostras

        FILE* pReadFile = fopen("benchmark.bin", "rb");
        if (!pReadFile) {
            printf("Erro: as amostras nao podem ser lidas.\n");
            return (2);
        }

        for (long i = 0; i < sample_size; i++) {
            // Leitura das amostras do arquivo
            if (i < raw_sample_size) {
            fread(&buffer, sizeof(double), 1, pReadFile);
            // Casting buffer as double (provavel de funcionar apenas se o .bin gerado vier da mesma máquina)
            converted_data = *((double*)buffer);
            samples[i] = converted_data;
            // Debug das amostras
            // printf("Var[%i]: %f\n", i, samples[i]); 
            } else { 
                // Preenchendo as posições com 0 caso o n de amostras não seja multiplo de 2
                samples[i] = 0;
            }
        }
        fclose(pReadFile);


        // ======================================================================
        //                                 FFT
        // ======================================================================
        start = clock();
        
        // ----------------------------Bit-Reversal------------------------------
        // Variáveis auxiliares (K e L serão utilizados mais a frente também)
        long K = 0, L = 1;
        complex aux = 0;
        // Loop para comparar todos os pares possíveis (N-1 pares)
        for (long i = 1; i < sample_size; i++) {
            if (i < L) { // Efetua a troca
                aux = samples[(L-1)];
                samples[(L-1)] = samples[(i-1)];
                samples[(i-1)] = aux;
            }
            K = sample_size / 2;
            while (K < L) {
                L -= K;
                K /= 2;
            }
            L += K;
        }
        
        // ------------------------Butterfly Operations--------------------------
        // Representações dos Twiddle factors (TF's) - Operadores de rotação Wn:
        // (U = Exponencial de varredura; W = Passo da varredura)
        complex U = 0 + 0*I, W = 0 + 0*I; 
        // Twiddle factor do butterfly (Wn a ser considerado no cálculo)
        complex T = 0 + 0*I;
        // L_exp = N*T (amostras*intervalo) --> Pontos na sub-DFT
        long L_exp = 0;
        // Ordem das amostras (em base 2)
        double M = log2(sample_size);
        // Índice ímpar
        long IQ = 0; 
        // IP = Índice par
        long IP = 0;

        // Fazendo M estágios de butterflies
        for (L = 1; L <= M; L++) {  
            L_exp = 1 << L; // L_exp = 2^L
            U = 1 + 0*I;
            // W = cexp((-2.0*M_PI) / L_exp); // Euler's representation
            W = cos((2.0*M_PI) / L_exp) - I*sin((2.0*M_PI) / L_exp); // Cos/Sin representation

            // Butterflies do estágio L-ésimo
            for (K = 0; K < (L_exp/2); K++) {

                // Computando butterflies que possuem simetria (usam o mesmo Wn)
                for (IP = K; IP < sample_size; IP += L_exp) {
                    // Indexação do vetor ímpar correspondente no butterfly
                    IQ = IP + (L_exp / 2);
                    // Cálculo do Twiddle Factor Wn 
                    T = samples[IQ] * U;
                    // Butterfly multiplication (explorando o princípio da simetria dos TF's)
                    samples[IQ] = samples[IP] - T;
                    samples[IP] = samples[IP] + T;
                }
                // Passo dos twiddle factors (Wn)
                U *= W;
            }
        }
        // O vetor samples[] agora tem os pontos calculados da fft.
        end = clock(); // Fim da FFT
        // Acumulando o número de ciclos de cada iteração:
        timer += (end - start);
    }
    
    // Profiling
    timer /= (double)CLOCKS_PER_SEC; // Tempo computado de 100 operações.
    printf("%f\n", timer);
    return 0;
}