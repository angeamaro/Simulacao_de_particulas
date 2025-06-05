#ifndef PARTICULAS_H
#define PARTICULAS_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1
#ifndef M_PI
#define M_PI 3.141592653589793
#endif

// Definição da estrutura da partícula
typedef struct {
    double x, y;      // Posição da partícula
    double vx, vy;    // Velocidade
    double massa;     // Massa
} Particula;

// Definição da estrutura da célula
typedef struct {
    double cx, cy;   // Centro de massa da célula
    double massa_total;  // Massa total das partículas dentro da célula
} Celula;

// Funções
void inicializar_particulas(int semente, double tamanho_espaco, int tamanho_grade, long long num_particulas, Particula *particulas);
void calcular_centros_de_massa(int grid_size, long long n_part, double side, Particula *par, Celula **grid);
long long calcular_forca_gravitacional(long long n_part, Particula *particulas);
long long simular_particulas(Particula *particulas, long long num_particulas, int num_passos, double tamanho_espaco, int grid_size, Celula **grid);
void print_result();
#endif // PARTICULAS_H