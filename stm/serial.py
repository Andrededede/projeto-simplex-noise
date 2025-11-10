import serial
import time
import threading
import sys
import numpy as np
# O 'cv2' não é mais necessário, removemos a dependência

# --- Configurações ---
SERIAL_PORT = 'COM7'  # Mude para a sua porta (ex: /dev/ttyUSB0 no Linux)
BAUD_RATE = 115200
IMG_WIDTH = 90
IMG_HEIGHT = 90
OUTPUT_BLOCK_LINES = 15 # 1/3 da imagem (30 linhas)
NUM_BLOCKS = 6

# Variável global para sinalizar o fim da thread
stop_sending = threading.Event()
raw_image_data = []

# --- Thread 1: O Transmissor (Sem modificações) ---
def sender_thread(ser):
    global raw_image_data
    print("[TX] Thread Transmissora iniciada.")
    
    # Parseia o PGM para uma lista de 90 bytearrays
    raw_image_data = []
    try:
        with open("entrada.pgm", 'r') as f:
            assert f.readline().strip() == 'P2' # Magia
            line = f.readline()
            while line.startswith('#'): line = f.readline() # Comentários
            w, h = map(int, line.split()) # Dimensões
            max_val = int(f.readline().strip()) # Max val
            pixels_str = f.read().split()
            pixels_int = [int(p) for p in pixels_str if p]
            
        for i in range(IMG_HEIGHT):
            linha = bytearray(pixels_int[i*IMG_WIDTH : (i+1)*IMG_WIDTH])
            raw_image_data.append(linha)
            
        print(f"[TX] Imagem de entrada parseada. {len(raw_image_data)} linhas.")
    except Exception as e:
        print(f"[TX] Erro ao ler 'entrada.pgm': {e}")
        print("[TX] Verifique se o arquivo está na pasta correta.")
        stop_sending.set() # Sinaliza para a thread main parar
        return

    line_num = 0
    while not stop_sending.is_set():
        try:
            # Cria o pacote: [NUM_LINHA (1 byte)][DADOS (90 bytes)]
            packet = bytes([line_num]) + raw_image_data[line_num]
            ser.write(packet) # Envia o pacote de 91 bytes

            sys.stdout.write(f"\rEnviando Linha: {line_num:02d}...")
            sys.stdout.flush()
            
            line_num = (line_num + 1) % IMG_HEIGHT # Loop de 0 a 89
            time.sleep(0.01) # Pausa
            
        except Exception as e:
            # Ocorre se a porta for fechada pela thread principal
            break
    print("[TX] Thread Transmissora finalizada.")


# --- NOVA FUNÇÃO ---
def save_as_pgm(filename, matrix):
    """
    Salva uma matriz NumPy (uint8) como um arquivo PGM P2 (ASCII).
    """
    print(f"\nSalvando matriz como '{filename}' no formato PGM P2...")
    try:
        # Pega as dimensões da matriz
        height, width = matrix.shape
        max_val = 255 # Valor máximo de pixel

        with open(filename, 'w') as f:
            # 1. Escreve o cabeçalho PGM
            f.write("P2\n")
            f.write(f"{width} {height}\n")
            f.write(f"{max_val}\n")

            # 2. Escreve os dados de pixel
            for y in range(height):
                # Pega todos os pixels da linha
                row_data = matrix[y]
                # Converte todos os números em strings
                row_strings = [str(pixel) for pixel in row_data]
                # Junta tudo com espaços
                f.write(" ".join(row_strings) + "\n")
        
        print(f"Arquivo '{filename}' salvo com sucesso.")

    except Exception as e:
        print(f"Erro ao salvar o arquivo PGM: {e}")


# --- Thread 2: O Receptor (Loop Principal Modificado) ---
def main():
    try:
        ser = serial.Serial(SERIAL_PORT, BAUD_RATE, timeout=5)
    except Exception as e:
        print(f"Erro ao abrir serial: {e}")
        return

    # A "matriz de inteiros" para salvar a saída
    output_matrix = np.zeros((IMG_HEIGHT, IMG_WIDTH), dtype=np.uint8)
    
    # Inicia a thread transmissora
    sender = threading.Thread(target=sender_thread, args=(ser,))
    sender.start()

    print("[RX] Aguardando dados do STM32 (printf)...")

    try:
        current_block = 0
        current_line_in_block = 0
        is_in_pgm_data = False

        while current_block < NUM_BLOCKS:
            # Se a thread de envio falhar (ex: não achou entrada.pgm), paramos.
            if stop_sending.is_set():
                print("[RX] Thread de envio falhou. Encerrando receptor.")
                break
                
            # 1. Lê uma linha de texto do terminal do STM32
            line = ser.readline().decode('ascii', errors='ignore').strip()

            if not line:
                # Timeout, provavelmente o STM32 está processando
                print("[RX] Aguardando STM32 processar...")
                continue

            # 2. Ignora mensagens de debug
            if line.startswith("---") or line.startswith("Processando Bloco"):
                print(f"STM32: {line}")
                continue
                
            # 3. Detecta o cabeçalho do PGM (início de um bloco)
            if line == "P2":
                print(f"[RX] Detectado início do Bloco {current_block + 1}")
                is_in_pgm_data = True
                current_line_in_block = 0
                # Lê e ignora as próximas duas linhas do cabeçalho (dimensões e max_val)
                ser.readline() 
                ser.readline()
                continue
                
            # 4. Se estamos dentro de um bloco, parseia os dados
            if is_in_pgm_data:
                try:
                    # Pega os números da linha (ex: "106 99 97...")
                    pixels = list(map(int, line.split()))
                    
                    if len(pixels) == IMG_WIDTH:
                        # Salva na nossa matriz de inteiros
                        y_global = (current_block * OUTPUT_BLOCK_LINES) + current_line_in_block
                        output_matrix[y_global, :] = pixels
                        
                        sys.stdout.write(f"\r[RX] Recebida linha de saída {y_global}...")
                        sys.stdout.flush()
                        
                        current_line_in_block += 1
                        
                        # Verifica se o bloco terminou
                        if current_line_in_block == OUTPUT_BLOCK_LINES:
                            is_in_pgm_data = False
                            current_block += 1
                            
                except Exception as e:
                    print(f"\nErro ao parsear linha de dados: {line} ({e})")
        
        if not stop_sending.is_set():
            print("\n[RX] Todos os blocos recebidos!")

    except KeyboardInterrupt:
        print("\n[RX] Interrompido pelo usuário.")
    except Exception as e:
        print(f"\n[RX] Erro na thread receptora: {e}")
        
    finally:
        # Sinaliza para a thread transmissora parar
        stop_sending.set()
        sender.join() # Espera a thread fechar
        ser.close()
        print("Conexão serial fechada.")

    # --- SEÇÃO DE SALVAMENTO (MODIFICADA) ---
    # Chama a nova função para salvar a matriz como PGM
    save_as_pgm("saida_distorcida.pgm", output_matrix)


if __name__ == "__main__":
    main()
