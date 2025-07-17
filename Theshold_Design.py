import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal

# ==== PARAMETERS ====
bit_rate = 1_000_000
fs = 42_000_000
samples_per_bit = fs // bit_rate
num_bits = 100000 # ใช้ค่ามากขึ้นเพื่อลดความผันผวนของ BER
epsilon = 1e-12
np.random.seed(42)

SNR_dB_range = np.arange(0, 11, .3)  # SNR = 0 ถึง 10 dB
alpha_values = [0.1,0.2,0.3, 0.5, 0.7, 0.9]  # EMA coefficient ที่ต้องการเปรียบเทียบ

ber_results = {}

# ==== LOOP สำหรับแต่ละค่า α ====
for alpha in alpha_values:
    ber_list = []
    for SNR_dB in SNR_dB_range:
        # สร้างข้อมูลบิตและ OOK signal
        data_bits = np.random.randint(0, 2, num_bits)
        ook_signal = np.repeat(data_bits, samples_per_bit).astype(float)

        # เพิ่ม noise แบบ AWGN
        signal_power = np.mean(ook_signal ** 2)
        SNR_linear = 10**(SNR_dB / 10)
        noise_power = signal_power / SNR_linear
        noise = normal(0, np.sqrt(noise_power), len(ook_signal))
        received_signal = ook_signal + noise

        # === LOG DETECTOR ===
        log_power = []
        for i in range(num_bits):
            segment = received_signal[i * samples_per_bit : (i + 1) * samples_per_bit]
            energy = np.sum(segment ** 2)
            log_energy = np.log10(energy + epsilon)
            log_power.append(log_energy)
        log_power = np.array(log_power)

        # === EMA ADAPTIVE THRESHOLD ===
        threshold = np.zeros(num_bits)
        threshold[0] = log_power[0]
        for i in range(1, num_bits):
            threshold[i] = alpha * log_power[i] + (1 - alpha) * threshold[i - 1]

        # === DECISION ===
        detected_bits = (log_power > threshold).astype(int)
        ber = np.sum(detected_bits != data_bits) / num_bits
        ber_list.append(ber)

    ber_results[f"α = {alpha}"] = ber_list

# ==== PLOT RESULT ====
plt.figure(figsize=(10, 6))
for label, ber_list in ber_results.items():
    plt.plot(SNR_dB_range, ber_list, marker='o', label=label)

plt.title("BER vs SNR for Different EMA α (Log Detector OOK)")
plt.xlabel("SNR (dB)")
plt.ylabel("Bit Error Rate (BER)")
plt.grid(True)
plt.legend()
plt.yscale("log")
plt.tight_layout()
plt.show()
