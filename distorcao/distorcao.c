#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ================================================================== */
/* IMPLEMENTAÇÃO 2D DO SIMPLEX NOISE (SEM MALLOC)            */
/* Baseado na implementação de domínio público de Stefan Gustavson      */
/* ================================================================== */

// --- Constantes estáticas (sem alocação) ---

// Tabela de permutação duplicada (512 elementos) para evitar % 256
static const unsigned char perm[256] = {
    151,160,137, 91, 90, 15,131, 13,201, 95, 96, 53,194,233, 7,225,
    140, 36,103, 30, 69,142, 8, 99, 37,240, 21, 10, 23,190, 6,148,
    247,120,234, 75, 0, 26,197, 62, 94,252,219,203,117, 35, 11, 32,
     57,177, 33, 88,237,149, 56, 87,174, 20,125,136,171,168, 68,175,
     74,165, 71,134,139, 48, 27,166, 77,146,158,231, 83,111,229,122,
     60,211,133,230,220,105, 92, 41, 55, 46,245, 40,244,102,143, 54,
     65, 25, 63,161, 1,216, 80, 73,209, 76,132,187,208, 89, 18,169,
    200,196,135,130,116,188,159, 86,164,100,109,198,173,186, 3, 64,
     52,217,226,250,124,123, 5,202, 38,147,118,126,255, 82, 85,212,
    207,206, 59,227, 47, 16, 58, 17,182,189, 28, 42,223,183,170,213,
    119,248,152, 2, 44,154,163, 70,221,153,101,155,167, 43,172, 9,
    129, 22, 39,253, 19, 98,108,110, 79,113,224,232,178,185,112,104,
    218,246, 97,228,251, 34,242,193,238,210,144, 12,191,179,162,241,
     81, 51,145,235,249, 14,239,107, 49,192,214, 31,181,199,106,157,
    184, 84,204,176,115,121, 50, 45,127, 4,150,254,138,236,205, 93,
    222,114, 67, 29, 24, 72,243,141,128,195, 78, 66,215, 61,156,180
};

// Tabela de gradientes 2D
static const int grad2[12][2] = {
    {1,1}, {-1,1}, {1,-1}, {-1,-1},
    {1,0}, {-1,0}, {1,0}, {-1,0}, // {1,0} e {-1,0} estão duplicados
    {0,1}, {0,-1}, {0,1}, {0,-1}  // {0,1} e {0,-1} estão duplicados
};

// --- Funções Auxiliares (inline) ---

// Floor rápido para ponto flutuante
static inline int fast_floor(float x) {
    return (int)floorf(x);
}

// Produto escalar 2D
static inline float dot(const int* g, float x, float y) {
    return g[0] * x + g[1] * y;
}

// --- Função Principal do Ruído 2D ---

float simplex_noise2d(float x_in, float y_in) {
    float n0, n1, n2; // Contribuições dos 3 cantos do simplexo

    // Fator de Skew (inclinação) da grade quadrada para a triangular
    // F2 = 0.5 * (sqrt(3.0) - 1.0)
    const float F2 = 0.366025403f;
    float s = (x_in + y_in) * F2; // Skew
    int i = fast_floor(x_in + s);
    int j = fast_floor(y_in + s);

    // Fator de Unskew (retorno) da grade triangular para a quadrada
    // G2 = (3.0 - sqrt(3.0)) / 6.0
    const float G2 = 0.211324865f;
    float t = (i + j) * G2;
    float X0 = i - t; // Coordenada X (unskewed) do canto da célula
    float Y0 = j - t; // Coordenada Y (unskewed) do canto da célula
    float x0 = x_in - X0; // Posição do ponto dentro da célula
    float y0 = y_in - Y0;

    // Determina qual dos dois triângulos da célula o ponto está
    int i1, j1; // Offsets para o segundo canto
    if (x0 > y0) {
        i1 = 1; j1 = 0; // Canto (i+1, j)
    } else {
        i1 = 0; j1 = 1; // Canto (i, j+1)
    }

    // Coordenadas relativas dos outros dois cantos
    float x1 = x0 - i1 + G2;
    float y1 = y0 - j1 + G2;
    float x2 = x0 - 1.0f + 2.0f * G2;
    float y2 = y0 - 1.0f + 2.0f * G2;

    // "Hash" das coordenadas dos cantos para obter índices da tabela de gradiente
    int ii = i % 256;
    int jj = j % 256;
    int gi0 = perm[(ii + perm[jj]) % 256] % 12;
    int gi1 = perm[(ii + i1 + perm[jj + j1] % 256)] % 12;
    int gi2 = perm[(ii + 1 + perm[jj + 1]) % 256] % 12;

    // Calcula a contribuição de cada canto
    float t0 = 0.5f - x0*x0 - y0*y0;
    if (t0 < 0.0f) {
        n0 = 0.0f;
    } else {
        t0 *= t0;
        n0 = t0 * t0 * dot(grad2[gi0], x0, y0);
    }

    float t1 = 0.5f - x1*x1 - y1*y1;
    if (t1 < 0.0f) {
        n1 = 0.0f;
    } else {
        t1 *= t1;
        n1 = t1 * t1 * dot(grad2[gi1], x1, y1);
    }

    float t2 = 0.5f - x2*x2 - y2*y2;
    if (t2 < 0.0f) {
        n2 = 0.0f;
    } else {
        t2 *= t2;
        n2 = t2 * t2 * dot(grad2[gi2], x2, y2);
    }

    // Soma as contribuições e escala o resultado para o intervalo [-1, 1]
    // (O fator 70.0 é empírico para normalizar)
    return 70.0f * (n0 + n1 + n2);
}

/* ================================================================== */
/* FIM DA IMPLEMENTAÇÃO DO SIMPLEX NOISE                 */
/* ================================================================== */

// --- Definições da Imagem ---
#define LARGURA 90
#define ALTURA 90
#define MAX_VALOR 255

// --- REQUISITO DO TRABALHO (8KB) ---
// A estrutura de dados de ENTRADA (aprox. 8KB)
unsigned char imagem_in[ALTURA][LARGURA];
// A estrutura de dados de SAÍDA
unsigned char imagem_out[ALTURA][LARGURA];


/* ===================================================================
 * PARTE 1: FUNÇÕES DE LEITURA E ESCRITA PGM P2
 * ===================================================================
 */

// Função para LER uma imagem PGM P2
void ler_imagem_pgm(const char* nome_arquivo) {
    FILE* fp;
    char magic[3];
    int largura_lida, altura_lida, max_lido;
    int pixel_val;

    fp = fopen(nome_arquivo, "r");
    if (fp == NULL) {
        perror("Erro ao abrir arquivo de entrada");
        fprintf(stderr, "Verifique se o arquivo '%s' existe na mesma pasta.\n", nome_arquivo);
        exit(1);
    }

    // 1. Ler cabeçalho
    fscanf(fp, "%2s", magic);
    if (strcmp(magic, "P2") != 0) {
        fprintf(stderr, "Erro: O arquivo de entrada nao e um PGM P2 (ASCII).\n");
        exit(1);
    }
    
    // Ignorar comentários (linhas que começam com #)
    char c = getc(fp);
    while (c == '\n' || c == '\r' || c == ' ' || c == '\t') c = getc(fp);
    if (c == '#') {
        while (c != '\n' && c != '\r') c = getc(fp);
    } else {
        ungetc(c, fp);
    }

    fscanf(fp, "%d %d", &largura_lida, &altura_lida);
    fscanf(fp, "%d", &max_lido);

    if (largura_lida != LARGURA || altura_lida != ALTURA) {
        fprintf(stderr, "Erro: A imagem de entrada tem tamanho %dx%d, mas o programa espera %dx%d.\n",
                largura_lida, altura_lida, LARGURA, ALTURA);
        exit(1);
    }

    // 2. Ler os dados de pixel
    for (int y = 0; y < ALTURA; y++) {
        for (int x = 0; x < LARGURA; x++) {
            if (fscanf(fp, "%d", &pixel_val) != 1) {
                 fprintf(stderr, "Erro ao ler os pixels da imagem.\n");
                 exit(1);
            }
            imagem_in[y][x] = (unsigned char)pixel_val;
        }
    }

    fclose(fp);
    printf("Imagem de entrada '%s' lida com sucesso.\n", nome_arquivo);
}

// Função para ESCREVER uma imagem PGM P2
void escrever_imagem_pgm(const char* nome_arquivo) {
    FILE* fp;
    fp = fopen(nome_arquivo, "w");
    if (fp == NULL) {
        perror("Erro ao criar arquivo de saida");
        exit(1);
    }

    // 1. Escrever o cabeçalho PGM P2
    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", LARGURA, ALTURA);
    fprintf(fp, "%d\n", MAX_VALOR);

    // 2. Escrever os dados de pixel
    for (int y = 0; y < ALTURA; y++) {
        for (int x = 0; x < LARGURA; x++) {
            fprintf(fp, "%d ", imagem_out[y][x]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    printf("Imagem de saida '%s' escrita com sucesso.\n", nome_arquivo);
}

/* ===================================================================
 * PARTE 2: O ALGORITMO DE DISTORÇÃO (OBJETIVO FINAL)
 * ===================================================================
 */
void processar_distorcao() {
    
    // Parâmetros para controlar o efeito
    float frequencia = 0.01f; // "Zoom" do ruído (quanto menor, mais suave)
    float intensidade = 10.0f; // "Força" da distorção (em pixels)

    printf("Processando distorcao com Simplex Noise...\n");

    // --- REQUISITO: Laços aninhados ---
    // Percorre cada pixel da imagem de SAÍDA
    for (int y = 0; y < ALTURA; y++) {
        for (int x = 0; x < LARGURA; x++) {
            
            // 1. Calcular o deslocamento (offset) usando o Simplex Noise
            // Gera um valor de ruído para o deslocamento em X
            float offsetX = simplex_noise2d(x * frequencia, y * frequencia);
            
            // Gera um segundo valor de ruído (em um "plano" 3D diferente) para o Y
            // Isso garante que o deslocamento X e Y não sejam correlacionados
            float offsetY = simplex_noise2d((x + 52.3f) * frequencia, (y - 13.7f) * frequencia);

            // 2. Calcular as novas coordenadas de amostragem
            // Para onde esse pixel (x,y) de saída deve "olhar" na entrada?
            int newX = (int)floorf(x + (offsetX * intensidade));
            int newY = (int)floorf(y + (offsetY * intensidade));

            // 3. Validação de Limites (Clamp)
            // Garante que não vamos tentar ler fora da memória da imagem_in
            if (newX < 0) newX = 0;
            if (newX >= LARGURA) newX = LARGURA - 1;
            if (newY < 0) newY = 0;
            if (newY >= ALTURA) newY = ALTURA - 1;

            // 4. Amostrar da imagem de entrada e escrever na saída
            // O pixel (x,y) da saída recebe a cor do pixel (newX, newY) da entrada.
            imagem_out[y][x] = imagem_in[newY][newX];
        }
    }
}

/* ===================================================================
 * PARTE 3: FUNÇÃO PRINCIPAL (Validação no PC)
 * ===================================================================
 */
int main() {
    
    // 1. Ler a imagem de entrada (deve se chamar "entrada.pgm")
    ler_imagem_pgm("entrada.pgm");

    // 2. Executar o algoritmo de distorção
    processar_distorcao();

    // 3. Escrever o resultado no arquivo PGM
    escrever_imagem_pgm("saida_distorcida.pgm");

    printf("Processo concluido.\n");
    return 0;
}