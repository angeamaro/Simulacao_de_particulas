#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "particulas.h"

int main(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Uso: ./parsim <semente> <tamanho do espaço> <tamanho da grade> <número de partículas> <número de passos>\n");
        return 1;
    }

    // Conversão dos argumentos
    int semente = atoi(argv[1]);
    double tamanho_espaco = atof(argv[2]);
    int tamanho_grade = atoi(argv[3]);  // Tamanho da grade
    long long num_particulas = atoll(argv[4]);
    int num_passos = atoi(argv[5]);

    // Verificação de valores válidos
    if (semente == 0 || tamanho_espaco <= 0 || tamanho_grade < 3 || num_particulas <= 0 || num_passos <= 0) {
        fprintf(stderr, "Erro: Parâmetros inválidos. Certifique-se de fornecer valores positivos válidos.\n");
        return 1;
    }

    // Alocar memória para as partículas
    Particula *particulas = (Particula *)malloc(num_particulas * sizeof(Particula));
    if (!particulas) {
        fprintf(stderr, "Erro ao alocar memória para partículas.\n");
        return 1;
    }

    // Alocar memória para a grade de células
    Celula **grade = (Celula **)malloc(tamanho_grade * sizeof(Celula *));
    if (!grade) {
        fprintf(stderr, "Erro ao alocar memória para a grade de células.\n");
        free(particulas);
        return 1;
    }

    for (int i = 0; i < tamanho_grade; i++) {
        grade[i] = (Celula *)malloc(tamanho_grade * sizeof(Celula));
        if (!grade[i]) {
            fprintf(stderr, "Erro ao alocar memória para a linha %d da grade.\n", i);
            for (int j = 0; j < i; j++) {
                free(grade[j]);
            }
            free(grade);
            free(particulas);
            return 1;
        }
    }

    // Inicializar partículas
    inicializar_particulas(semente, tamanho_espaco, tamanho_grade, num_particulas, particulas);
      double tempo = -omp_get_wtime();
    // Simular as partículas, passando todos os parâmetros necessários
    simular_particulas(particulas, num_particulas, num_passos, tamanho_espaco, tamanho_grade, grade);
    tempo += omp_get_wtime();
    printf("Tempo total: %f segundos\n", tempo);
    // Liberar memória alocada
    for (int i = 0; i < tamanho_grade; i++) {
        free(grade[i]);
    }
    free(grade);
    free(particulas);

    return 0;
}
