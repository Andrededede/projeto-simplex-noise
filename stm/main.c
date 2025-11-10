/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.c
  * @brief          : Main program body
  ******************************************************************************
  */
/* USER CODE END Header */
/* Includes ------------------------------------------------------------------*/
#include "main.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

#include <stdio.h>  // Para printf
#include <math.h>   // Para floorf
#include <string.h> // Para memset e memcpy

/* USER CODE END Includes */

/* Private typedef -----------------------------------------------------------*/
/* USER CODE BEGIN PTD */

/* USER CODE END PTD */

/* Private define ------------------------------------------------------------*/
/* USER CODE BEGIN PD */

// --- Configurações da Arquitetura (CORRIGIDO) ---
#define IMG_WIDTH 90
#define IMG_HEIGHT 90
#define MAX_VALOR 255

// Cache de Entrada: 1/2 da imagem (45 linhas)
#define CACHE_IN_LINES 45
// Buffer de Saída: 15 linhas
#define OUTPUT_BLOCK_LINES 15 // CORRIGIDO
// Número de blocos de saída
#define NUM_OUTPUT_BLOCKS 6 // CORRIGIDO (90 / 15 = 6)

/* USER CODE END PD */

/* Private macro -------------------------------------------------------------*/
/* USER CODE BEGIN PM */

/* USER CODE END PM */

/* Private variables ---------------------------------------------------------*/
UART_HandleTypeDef huart2;

/* USER CODE BEGIN PV */

// --- Buffers Globais (na RAM de 8KB) ---
uint8_t cache_in_data[CACHE_IN_LINES][IMG_WIDTH]; // 45 * 90 = 4050 bytes
uint8_t cache_in_tags[CACHE_IN_LINES];            // 45 bytes (para o rótulo da linha 0-89)

// Buffer de saída agora tem 15 linhas
uint8_t output_block[OUTPUT_BLOCK_LINES][IMG_WIDTH]; // 15 * 90 = 1350 bytes

// --- Variáveis de Controle ---
volatile uint8_t g_uart_rx_buffer[IMG_WIDTH + 1];
volatile int g_cache_replace_index = 0;

/* USER CODE END PV */

/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_USART2_UART_Init(void);
/* USER CODE BEGIN PFP */

/* USER CODE END PFP */

/* Private user code ---------------------------------------------------------*/
/* USER CODE BEGIN 0 */

// --- Simplex Noise (Idêntico ao seu código original) ---
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
static const int grad2[12][2] = {
    {1,1}, {-1,1}, {1,-1}, {-1,-1},
    {1,0}, {-1,0}, {1,0}, {-1,0},
    {0,1}, {0,-1}, {0,1}, {0,-1}
};

// --- Funções Auxiliares (inline) ---
static inline int fast_floor(float x) {
    return (int)floorf(x);
}
static inline float dot(const int* g, float x, float y) {
    return g[0] * x + g[1] * y;
}

// --- Função Simplex Noise (Idêntica à sua) ---
float simplex_noise2d(float x_in, float y_in) {
    float n0, n1, n2;
    const float F2 = 0.366025403f;
    float s = (x_in + y_in) * F2;
    int i = fast_floor(x_in + s);
    int j = fast_floor(y_in + s);
    const float G2 = 0.211324865f;
    float t = (i + j) * G2;
    float X0 = i - t;
    float Y0 = j - t;
    float x0 = x_in - X0;
    float y0 = y_in - Y0;
    int i1, j1;
    if (x0 > y0) { i1 = 1; j1 = 0; }
    else { i1 = 0; j1 = 1; }
    float x1 = x0 - i1 + G2;
    float y1 = y0 - j1 + G2;
    float x2 = x0 - 1.0f + 2.0f * G2;
    float y2 = y0 - 1.0f + 2.0f * G2;
    int ii = i % 256;
    int jj = j % 256;
    int gi0 = perm[(ii + perm[jj]) % 256] % 12;
    int gi1 = perm[(ii + i1 + perm[jj + j1] % 256)] % 12;
    int gi2 = perm[(ii + 1 + perm[jj + 1]) % 256] % 12;
    float t0 = 0.5f - x0*x0 - y0*y0;
    if (t0 < 0.0f) { n0 = 0.0f; }
    else { t0 *= t0; n0 = t0 * t0 * dot(grad2[gi0], x0, y0); }
    float t1 = 0.5f - x1*x1 - y1*y1;
    if (t1 < 0.0f) { n1 = 0.0f; }
    else { t1 *= t1; n1 = t1 * t1 * dot(grad2[gi1], x1, y1); }
    float t2 = 0.5f - x2*x2 - y2*y2;
    if (t2 < 0.0f) { n2 = 0.0f; }
    else { t2 *= t2; n2 = t2 * t2 * dot(grad2[gi2], x2, y2); }
    return 70.0f * (n0 + n1 + n2);
}


// --- LÓGICA DE INTERRUPÇÃO E CACHE (ASSÍNCRONO) ---

/**
 * @brief Esta é a Interrupção! Chamada pela HAL quando 91 bytes chegam.
 */
void HAL_UART_RxCpltCallback(UART_HandleTypeDef *huart) {
    if (huart->Instance == USART2) { // Verifique se é a sua UART

        int received_line_num = g_uart_rx_buffer[0];
        uint8_t* data = (uint8_t*)&g_uart_rx_buffer[1];

        memcpy(cache_in_data[g_cache_replace_index], data, IMG_WIDTH);
        cache_in_tags[g_cache_replace_index] = received_line_num;

        g_cache_replace_index = (g_cache_replace_index + 1) % CACHE_IN_LINES;

        HAL_UART_Receive_IT(&huart2, (uint8_t*)g_uart_rx_buffer, IMG_WIDTH + 1);
    }
}

/**
 * @brief Verifica se uma linha (y_in) está no cache.
 * @return 1 se está no cache (Hit), 0 se não está (Miss).
 */
int is_line_in_cache(int y_in) {
    for (volatile int i = 0; i < CACHE_IN_LINES; i++) {
        if (cache_in_tags[i] == y_in) {
            return 1; // CACHE HIT!
        }
    }
    return 0; // CACHE MISS!
}

/**
 * @brief Pega o pixel do cache (DEPOIS de saber que ele existe).
 */
uint8_t get_pixel_from_cache(int y_in, int x_in) {
     for (int i = 0; i < CACHE_IN_LINES; i++) {
        if (cache_in_tags[i] == y_in) {
            return cache_in_data[i][x_in];
        }
    }
    return 0; // (Nunca deve chegar aqui)
}


// --- LÓGICA DE PROCESSAMENTO (COM BUSY-WAIT) ---

/**
 * @brief Processa UM bloco de saída (agora com 15 linhas)
 */
void processar_bloco_de_saida(int block_num) {

    // Agora é sempre 15 linhas
    int start_y = block_num * OUTPUT_BLOCK_LINES;
    int end_y = start_y + OUTPUT_BLOCK_LINES;

    float frequencia = 0.01f;
    float intensidade = 10.0f;

    for (int y = start_y; y < end_y; y++) {
        for (int x = 0; x < IMG_WIDTH; x++) {

            float offsetX = simplex_noise2d(x * frequencia, y * frequencia);
            float offsetY = simplex_noise2d((x + 52.3f) * frequencia, (y - 13.7f) * frequencia);
            int newX = (int)floorf(x + (offsetX * intensidade));
            int newY = (int)floorf(y + (offsetY * intensidade));

            if (newX < 0) newX = 0;
            if (newX >= IMG_WIDTH) newX = IMG_WIDTH - 1;
            if (newY < 0) newY = 0;
            if (newY >= IMG_HEIGHT) newY = IMG_HEIGHT - 1;

            // "Busy-Wait" (Espera Ocupada)
            while (is_line_in_cache(newY) == 0) {
                 HAL_Delay(1); // Pausa para a Interrupção trabalhar
            }

            uint8_t pixel = get_pixel_from_cache(newY, newX);

            // Salva no buffer de saída
            output_block[y - start_y][x] = pixel;
        }
    }
}

/**
 * @brief "Printa" o bloco de saída no terminal (formato PGM)
 */
void enviar_bloco_pelo_terminal(int block_num) {

    // Cabeçalho PGM (sempre 15 linhas)
    printf("P2\n");
    printf("%d %d\n", IMG_WIDTH, OUTPUT_BLOCK_LINES);
    printf("%d\n", MAX_VALOR);

    // Envia as 15 linhas
    for (int y = 0; y < OUTPUT_BLOCK_LINES; y++) {
        for (int x = 0; x < IMG_WIDTH; x++) {
            printf("%d ", output_block[y][x]);
        }
        printf("\n");
    }
    printf("\n--- FIM DO BLOCO %d ---\n\n", block_num + 1);
    memset(output_block, 0, sizeof(output_block));
}


// --- Função Principal (em main.c) ---
int main(void) {
    HAL_Init();
    SystemClock_Config();
    MX_GPIO_Init();
    MX_USART2_UART_Init();

    memset(cache_in_tags, 255, sizeof(cache_in_tags));

    printf("STM32 Pronto. Iniciando escuta assincrona...\n");

    // 1. Inicia a "escuta" da UART por Interrupção (modo Assíncrono)
    HAL_UART_Receive_IT(&huart2, (uint8_t*)g_uart_rx_buffer, IMG_WIDTH + 1);

    // 2. Loop principal (WHILE 1)
    while (1) {
        printf("Processando %d blocos...\n", NUM_OUTPUT_BLOCKS);
        HAL_Delay(10); // Pausa para o printf enviar

        // Loop principal (CORRIGIDO para 6 blocos)
        for (int i = 0; i < NUM_OUTPUT_BLOCKS; i++) {

            printf("--- Processando Bloco de Saida %d/%d... ---\n", i + 1, NUM_OUTPUT_BLOCKS);
            HAL_Delay(10);

            // 1. Processa (e espera por "Cache Miss" no busy-wait)
            processar_bloco_de_saida(i);

            // 2. Envia o bloco de saída pelo terminal
            enviar_bloco_pelo_terminal(i);
        }

        printf("\n\n--- PROCESSAMENTO CONCLUIDO --- \n");
        printf("Reiniciando em 10 segundos...\n\n");
        HAL_Delay(10000);
    }
}
/* USER CODE END 3 */

/**
  * @brief System Clock Configuration
  * @retval None
  */
void SystemClock_Config(void)
{
  RCC_OscInitTypeDef RCC_OscInitStruct = {0};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {0};

  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_ON;
  RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  RCC_OscInitStruct.PLL.PLLMUL = RCC_PLL_MUL12;
  RCC_OscInitStruct.PLL.PREDIV = RCC_PREDIV_DIV1;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK)
  {
    Error_Handler();
  }

  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK|RCC_CLOCKTYPE_SYSCLK
                              |RCC_CLOCKTYPE_PCLK1;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_1) != HAL_OK)
  {
    Error_Handler();
  }
}

/**
  * @brief USART2 Initialization Function
  * @param None
  * @retval None
  */
static void MX_USART2_UART_Init(void)
{

  /* USER CODE BEGIN USART2_Init 0 */

  /* USER CODE END USART2_Init 0 */

  /* USER CODE BEGIN USART2_Init 1 */

  /* USER CODE END USART2_Init 1 */
  huart2.Instance = USART2;
  huart2.Init.BaudRate = 115200;
  huart2.Init.WordLength = UART_WORDLENGTH_8B;
  huart2.Init.StopBits = UART_STOPBITS_1;
  huart2.Init.Parity = UART_PARITY_NONE;
  huart2.Init.Mode = UART_MODE_TX_RX;
  huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  huart2.Init.OverSampling = UART_OVERSAMPLING_16;
  huart2.Init.OneBitSampling = UART_ONE_BIT_SAMPLE_DISABLE;
  huart2.AdvancedInit.AdvFeatureInit = UART_ADVFEATURE_NO_INIT;
  if (HAL_UART_Init(&huart2) != HAL_OK)
  {
    Error_Handler();
  }
  /* USER CODE BEGIN USART2_Init 2 */

  /* USER CODE END USART2_Init 2 */

}

/**
  * @brief GPIO Initialization Function
  * @param None
  * @retval None
  */
static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};
  /* USER CODE BEGIN MX_GPIO_Init_1 */

  /* USER CODE END MX_GPIO_Init_1 */

  /* GPIO Ports Clock Enable */
  __HAL_RCC_GPIOC_CLK_ENABLE();
  __HAL_RCC_GPIOF_CLK_ENABLE();
  __HAL_RCC_GPIOA_CLK_ENABLE();

  /*Configure GPIO pin Output Level */
  HAL_GPIO_WritePin(LD2_GPIO_Port, LD2_Pin, GPIO_PIN_RESET);

  /*Configure GPIO pin : B1_Pin */
  GPIO_InitStruct.Pin = B1_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_IT_FALLING;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  HAL_GPIO_Init(B1_GPIO_Port, &GPIO_InitStruct);

  /*Configure GPIO pin : LD2_Pin */
  GPIO_InitStruct.Pin = LD2_Pin;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(LD2_GPIO_Port, &GPIO_InitStruct);

  /* USER CODE BEGIN MX_GPIO_Init_2 */

  /* USER CODE END MX_GPIO_Init_2 */
}

/* USER CODE BEGIN 4 */

/**
 * @brief Redireciona a função printf da biblioteca C para a UART.
 */
#ifdef __GNUC__
/* Com GCC, sobrescrevemos _write. */
int _write(int file, char *ptr, int len)
{
  (void)file; // Evita warning de 'unused'
  HAL_UART_Transmit(&huart2, (uint8_t*)ptr, len, HAL_MAX_DELAY); // Timeout longo
  return len;
}
#else
/* Com Keil/IAR, sobrescrevemos fputc */
int fputc(int ch, FILE *f)
{
  HAL_UART_Transmit(&huart2, (uint8_t *)&ch, 1, HAL_MAX_DELAY);
  return ch;
}
#endif

/* USER CODE END 4 */

/**
  * @brief  This function is executed in case of error occurrence.
  * @retval None
  */
void Error_Handler(void)
{
  /* USER CODE BEGIN Error_Handler_Debug */
  /* User can add his own implementation to report the HAL error return state */
  __disable_irq();
  while (1)
  {
  }
  /* USER CODE END Error_Handler_Debug */
}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  * where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     ex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif /* USE_FULL_ASSERT */
