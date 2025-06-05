#include "particulas.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned int semente_aleatoria;

// Inicializa o gerador com uma semente personalizada
void iniciar_gerador_aleatorio(int entrada) {
    semente_aleatoria = entrada + 987654321;
}

// Retorna um número aleatório uniformemente distribuído entre 0 e 1
double numero_aleatorio_uniforme() {
    int semente_temp = semente_aleatoria;
    semente_aleatoria ^= (semente_aleatoria << 13);
    semente_aleatoria ^= (semente_aleatoria >> 17);
    semente_aleatoria ^= (semente_aleatoria << 5);
    return 0.5 + 0.2328306e-09 * (semente_temp + (int)semente_aleatoria);
}

// Retorna um número aleatório com distribuição normal truncada em [0,1)
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

// Inicializa as partículas com posições, velocidades e massas
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

// Função auxiliar para garantir que valores estão no intervalo [0, lado)
static inline double normalizar_posicao(double pos, double lado) {
    pos = fmod(pos, lado);
    return pos < 0 ? pos + lado : pos;
}

// Função auxiliar para obter o índice da célula
static inline int indice_celula(double coord, double lado, int grid_size) {
    int idx = (int)(coord / (lado / grid_size));
    if (idx < 0) return 0;
    if (idx >= grid_size) return grid_size - 1;
    return idx;
}

// Calcula os centros de massa de cada célula da grade
void calcular_centros_de_massa(int grid_size, long long n_part, double lado, Particula *par, Celula **grid) {
    // Zera os valores da grade
    for (int i = 0; i < grid_size; i++)
        for (int j = 0; j < grid_size; j++)
            grid[i][j] = (Celula){0, 0, 0};

    // Acumula massa e posição ponderada por massa
    for (long long i = 0; i < n_part; i++) {
        if (par[i].massa <= 0) continue;

        double x = normalizar_posicao(par[i].x, lado);
        double y = normalizar_posicao(par[i].y, lado);
        int cx = indice_celula(x, lado, grid_size);
        int cy = indice_celula(y, lado, grid_size);

        Celula *cel = &grid[cx][cy];
        cel->massa_total += par[i].massa;
        cel->cx += par[i].x * par[i].massa;
        cel->cy += par[i].y * par[i].massa;
    }

    // Finaliza o cálculo dos centros de massa
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

// Calcula a força gravitacional entre as partículas
void calcular_forca_gravitacional(long long n_part, Particula *p, double lado, int grid_size, Celula **grid) {
    for (long long i = 0; i < n_part; i++) {
        if (p[i].massa <= 0) continue;

        double fx = 0, fy = 0;

        // Índice da célula onde a partícula i está
        int cx = indice_celula(p[i].x, lado, grid_size);
        int cy = indice_celula(p[i].y, lado, grid_size);

        // Interações com partículas da mesma célula
        for (long long j = 0; j < n_part; j++) {
            if (i == j || p[j].massa <= 0) continue;

            int cxj = indice_celula(p[j].x, lado, grid_size);
            int cyj = indice_celula(p[j].y, lado, grid_size);

            if (cx == cxj && cy == cyj) {
                double dx = p[j].x - p[i].x;
                double dy = p[j].y - p[i].y;

                // Ajustes toroidais
                if (dx > 0.5 * lado) dx -= lado;
                if (dx < -0.5 * lado) dx += lado;
                if (dy > 0.5 * lado) dy -= lado;
                if (dy < -0.5 * lado) dy += lado;

                double dist2 = dx * dx + dy * dy + EPSILON2;
                double inv_dist = 1.0 / sqrt(dist2);
                double forca = G * p[i].massa * p[j].massa / dist2;

                fx += forca * dx * inv_dist;
                fy += forca * dy * inv_dist;
            }
        }

        // Interações com os centros de massa das 8 células vizinhas
        for (int dx_cell = -1; dx_cell <= 1; dx_cell++) {
            for (int dy_cell = -1; dy_cell <= 1; dy_cell++) {
                if (dx_cell == 0 && dy_cell == 0) continue; // já tratamos a própria célula acima

                int nx = (cx + dx_cell + grid_size) % grid_size;
                int ny = (cy + dy_cell + grid_size) % grid_size;

                Celula *vizinha = &grid[nx][ny];
                if (vizinha->massa_total > 0) {
                    double dx = vizinha->cx - p[i].x;
                    double dy = vizinha->cy - p[i].y;

                    // Ajustes toroidais
                    if (dx > 0.5 * lado) dx -= lado;
                    if (dx < -0.5 * lado) dx += lado;
                    if (dy > 0.5 * lado) dy -= lado;
                    if (dy < -0.5 * lado) dy += lado;

                    double dist2 = dx * dx + dy * dy + EPSILON2;
                    double inv_dist = 1.0 / sqrt(dist2);
                    double forca = G * p[i].massa * vizinha->massa_total / dist2;

                    fx += forca * dx * inv_dist;
                    fy += forca * dy * inv_dist;
                }
            }
        }

        // Atualiza a velocidade
        p[i].vx += fx / p[i].massa * DELTAT;
        p[i].vy += fy / p[i].massa * DELTAT;
    }
}


// Simula o movimento das partículas por um número de passos
void simular_particulas(Particula *p, long long n, int passos, double lado, int grid_size, Celula **grid) {
    long long num_colisoes = 0;
    for (int t = 0; t < passos; t++) {
        calcular_centros_de_massa(grid_size, n, lado, p, grid);
        calcular_forca_gravitacional(n, p, lado, grid_size, grid);
;

        for (long long i = 0; i < n; i++) {
            if (p[i].massa <= 0) continue;

            // Atualiza posições com ajuste toroidal
            p[i].x = normalizar_posicao(p[i].x + p[i].vx * DELTAT, lado);
            p[i].y = normalizar_posicao(p[i].y + p[i].vy * DELTAT, lado);
        }

        // Detecção e fusão de colisões (sem loop aninhado redundante)
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
                    // Fusão de partículas
                    double m_total = p[i].massa + p[j].massa;
                    p[i].x = (p[i].x * p[i].massa + p[j].x * p[j].massa) / m_total;
                    p[i].y = (p[i].y * p[i].massa + p[j].y * p[j].massa) / m_total;
                    p[i].vx = (p[i].vx * p[i].massa + p[j].vx * p[j].massa) / m_total;
                    p[i].vy = (p[i].vy * p[i].massa + p[j].vy * p[j].massa) / m_total;
                    p[i].massa = m_total;

                    p[j].massa = 0; // Marca como fundida
                    num_colisoes++;
                }
            }
        }
    }

    // Apenas imprime a posição final da primeira partícula (como referência)
    printf("%.3f %.3f\n", p[0].x, p[0].y);

    printf("%lld\n", num_colisoes);


}
