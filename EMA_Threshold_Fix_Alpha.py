import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.random import normal

# ==== PARAMETERS ====
bit_rate = 1000_000                 # bits per second
fs = 100_000_000                     # sampling rate (Hz)
samples_per_bit = fs // bit_rate
num_bits = 50                  # number of bits
SNR_dB = 5                      # Signal-to-Noise Ratio in dB
alpha = 0.3            # EMA threshold coefficient
epsilon = 1e-12                 # to avoid log(0)
np.random.seed(42)             # for reproducibility

# ==== GENERATE RANDOM BITS ====
data_bits = np.random.randint(0, 2, num_bits)  # 0 or 1

# ==== OOK MODULATION ====
ook_signal = np.repeat(data_bits, samples_per_bit).astype(float)

# ==== ADD NOISE ====
signal_power = np.mean(ook_signal ** 2)
SNR_linear = 10**(SNR_dB / 10)
noise_power = signal_power / SNR_linear
noise = normal(0, np.sqrt(noise_power), len(ook_signal))
received_signal = ook_signal + noise

# ==== LOG DETECTOR ====
log_power = []
for i in range(num_bits):
    segment = received_signal[i * samples_per_bit : (i + 1) * samples_per_bit]
    energy = np.sum(segment ** 2)
    log_energy = np.log10(energy + epsilon)
    log_power.append(log_energy)
log_power = np.array(log_power)

# ==== ADAPTIVE THRESHOLD (EMA) ====
threshold = np.zeros(num_bits)
threshold[0] = log_power[0]
for i in range(1, num_bits):
    threshold[i] = alpha * log_power[i] + (1 - alpha) * threshold[i - 1]

# ==== DECISION ====
detected_bits = (log_power > threshold).astype(int)
ber = np.sum(detected_bits != data_bits) / num_bits

# ==== SAVE RESULTS ====
results = pd.DataFrame({
    "Original Bit": data_bits,
    "Log-Power": log_power,
    "Threshold": threshold,
    "Detected Bit": detected_bits
})
results.to_csv("log_detector_ook_results.csv", index=False)
print("📁 ผลลัพธ์บันทึกไว้ใน log_detector_ook_results.csv")
print(f"📊 Bit Error Rate (BER) = {ber:.2%}")

# ==== PLOT ====
plt.figure(figsize=(12, 6))
plt.plot(log_power, label="Log Power (Detector Output)", marker='o')
plt.plot(threshold, label="Adaptive Threshold (EMA)", linestyle='--')
plt.plot(data_bits, label="Original Bit", linestyle=':', alpha=0.5)
plt.plot(detected_bits, label="Detected Bit", linestyle='-.', alpha=0.5)
plt.title(f"OOK Detection with Log Detector + EMA Threshold\nSNR = {SNR_dB} dB, BER = {ber:.2%}")
plt.xlabel("Bit Index")
plt.ylabel("Log Power / Bit Value")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
