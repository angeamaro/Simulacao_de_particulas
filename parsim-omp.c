#include "particulas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define G 6.67408e-11
#define EPSILON2 (0.005 * 0.005)
#define DELTAT 0.1

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

unsigned int semente_aleatoria;

void iniciar_gerador_aleatorio(int entrada) {
    semente_aleatoria = entrada + 987654321;
}

double numero_aleatorio_uniforme() {
    int semente_temp = semente_aleatoria;
    semente_aleatoria ^= (semente_aleatoria << 13);
    semente_aleatoria ^= (semente_aleatoria >> 17);
    semente_aleatoria ^= (semente_aleatoria << 5);
    return 0.5 + 0.2328306e-09 * (semente_temp + (int)semente_aleatoria);
}

double numero_aleatorio_normal() {
    double u1, u2, z;
    do {
        u1 = numero_aleatorio_uniforme();
        u2 = numero_aleatorio_uniforme();
        z = sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        z = 0.5 + 0.15 * z;
    } while (z < 0 || z >= 1);
    return z;
}

void inicializar_particulas(int semente, double tamanho_espaco, int tamanho_grade, long long num_particulas, Particula *particulas) {
    double (*gerador)() = semente < 0 ? numero_aleatorio_normal : numero_aleatorio_uniforme;
    iniciar_gerador_aleatorio(abs(semente));

    double fator_espaco = tamanho_espaco / tamanho_grade / 5.0;
    double fator_massa = 0.01 * tamanho_grade * tamanho_grade / (num_particulas * G) * EPSILON2;

    for (long long i = 0; i < num_particulas; i++) {
        particulas[i] = (Particula){
            .x = gerador() * tamanho_espaco,
            .y = gerador() * tamanho_espaco,
            .vx = (gerador() - 0.5) * fator_espaco,
            .vy = (gerador() - 0.5) * fator_espaco,
            .massa = gerador() * fator_massa
        };
    }
}

static inline double normalizar_posicao(double pos, double lado) {
    pos = fmod(pos, lado);
    return pos < 0 ? pos + lado : pos;
}

static inline int indice_celula(double coord, double lado, int grid_size) {
    int idx = (int)(coord / (lado / grid_size));
    if (idx < 0) return 0;
    if (idx >= grid_size) return grid_size - 1;
    return idx;
}

void calcular_centros_de_massa(int grid_size, long long n_part, double lado, Particula *par, Celula **grid) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < grid_size; i++)
        for (int j = 0; j < grid_size; j++)
            grid[i][j] = (Celula){0, 0, 0};

    #pragma omp parallel for
    for (long long i = 0; i < n_part; i++) {
        if (par[i].massa <= 0) continue;

        double x = normalizar_posicao(par[i].x, lado);
        double y = normalizar_posicao(par[i].y, lado);
        int cx = indice_celula(x, lado, grid_size);
        int cy = indice_celula(y, lado, grid_size);

        #pragma omp atomic
        grid[cx][cy].massa_total += par[i].massa;
        #pragma omp atomic
        grid[cx][cy].cx += par[i].x * par[i].massa;
        #pragma omp atomic
        grid[cx][cy].cy += par[i].y * par[i].massa;
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            Celula *cel = &grid[i][j];
            if (cel->massa_total > 0) {
                cel->cx = normalizar_posicao(cel->cx / cel->massa_total, lado);
                cel->cy = normalizar_posicao(cel->cy / cel->massa_total, lado);
            }
        }
    }
}

long long calcular_forca_gravitacional(long long n_part, Particula *p) {
    #pragma omp parallel for
    for (long long i = 0; i < n_part; i++) {
        if (p[i].massa <= 0) continue;

        double fx = 0, fy = 0;

        for (long long j = 0; j < n_part; j++) {
            if (i == j || p[j].massa <= 0) continue;

            double dx = p[j].x - p[i].x;
            double dy = p[j].y - p[i].y;

            if (dx > 0.5) dx -= 1.0;
            if (dx < -0.5) dx += 1.0;
            if (dy > 0.5) dy -= 1.0;
            if (dy < -0.5) dy += 1.0;

            double dist2 = dx * dx + dy * dy + EPSILON2;
            double inv_dist = 1.0 / sqrt(dist2);
            double forca = G * p[i].massa * p[j].massa / dist2;

            fx += forca * dx * inv_dist;
            fy += forca * dy * inv_dist;
        }

        p[i].vx += fx / p[i].massa * DELTAT;
        p[i].vy += fy / p[i].massa * DELTAT;
    }

    return 0; // sem uso aqui, pois colisões são contadas depois
}

long long simular_particulas(Particula *p, long long n, int passos, double lado, int grid_size, Celula **grid) {
    long long num_colisoes = 0;

    for (int t = 0; t < passos; t++) {
        calcular_centros_de_massa(grid_size, n, lado, p, grid);
        calcular_forca_gravitacional(n, p);

        #pragma omp parallel for
        for (long long i = 0; i < n; i++) {
            if (p[i].massa <= 0) continue;
            p[i].x = normalizar_posicao(p[i].x + p[i].vx * DELTAT, lado);
            p[i].y = normalizar_posicao(p[i].y + p[i].vy * DELTAT, lado);
        }

        #pragma omp parallel for reduction(+:num_colisoes)
        for (long long i = 0; i < n; i++) {
            if (p[i].massa <= 0) continue;
            for (long long j = i + 1; j < n; j++) {
                if (p[j].massa <= 0) continue;

                double dx = p[j].x - p[i].x;
                double dy = p[j].y - p[i].y;

                if (dx > lado / 2) dx -= lado;
                if (dx < -lado / 2) dx += lado;
                if (dy > lado / 2) dy -= lado;
                if (dy < -lado / 2) dy += lado;

                if (dx * dx + dy * dy < EPSILON2) {
                    double m_total = p[i].massa + p[j].massa;
                    #pragma omp critical
                    {
                        if (p[i].massa > 0 && p[j].massa > 0) {
                            p[i].x = (p[i].x * p[i].massa + p[j].x * p[j].massa) / m_total;
                            p[i].y = (p[i].y * p[i].massa + p[j].y * p[j].massa) / m_total;
                            p[i].vx = (p[i].vx * p[i].massa + p[j].vx * p[j].massa) / m_total;
                            p[i].vy = (p[i].vy * p[i].massa + p[j].vy * p[j].massa) / m_total;
                            p[i].massa = m_total;
                            p[j].massa = 0;
                            num_colisoes++;
                        }
                    }
                }
            }
        }
    }

    return num_colisoes;
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Uso: %s <semente> <lado> <tamanho_grade> <n_particulas> <passos>\n", argv[0]);
        return 1;
    }

    int semente = atoi(argv[1]);
    double lado = atof(argv[2]);
    int grid_size = atoi(argv[3]);
    long long n_particulas = atoll(argv[4]);
    int passos = atoi(argv[5]);

    Particula *particulas = malloc(n_particulas * sizeof(Particula));

    Celula **grid = malloc(grid_size * sizeof(Celula *));
    for (int i = 0; i < grid_size; i++)
        grid[i] = malloc(grid_size * sizeof(Celula));

    inicializar_particulas(semente, lado, grid_size, n_particulas, particulas);

    double tempo = -omp_get_wtime();
    long long colisoes = simular_particulas(particulas, n_particulas, passos, lado, grid_size, grid);
    tempo += omp_get_wtime();

    // Saída padrão (obrigatória)
    printf("%.3f %.3f\n", particulas[0].x, particulas[0].y);
    printf("%lld\n", colisoes);

    // Tempo (stderr)
    fprintf(stderr, "%.1fs\n", tempo);

    // Liberação de memória
    for (int i = 0; i < grid_size; i++) free(grid[i]);
    free(grid);
    free(particulas);

    return 0;
}
