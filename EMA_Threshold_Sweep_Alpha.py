import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from numpy.random import normal

# ==== PARAMETERS ====
bit_rate = 5_000_000
fs = 100_000_000
samples_per_bit = fs // bit_rate
num_bits = 50
SNR_dB = 1
epsilon = 1e-12
np.random.seed(42)

data_bits = np.random.randint(0, 2, num_bits)
ook_signal = np.repeat(data_bits, samples_per_bit).astype(float)
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

# ==== TEST MULTIPLE ALPHA ====
alpha_list = [0.3, 0.5, 0.7, 0.9, 0.95]
ber_results = {}

plt.figure(figsize=(12, 6))

for alpha in alpha_list:
    threshold = np.zeros(num_bits)
    threshold[0] = log_power[0]
    for i in range(1, num_bits):
        threshold[i] = alpha * log_power[i] + (1 - alpha) * threshold[i - 1]

    detected_bits = (log_power > threshold).astype(int)
    ber = np.sum(detected_bits != data_bits) / num_bits
    ber_results[alpha] = ber

    plt.plot(threshold, label=f"Threshold α={alpha:.2f} (BER={ber:.2%})")

# ==== PLOT ====
plt.plot(log_power, 'k--', label="Log Power")
plt.plot(data_bits, ':', label="Original Bits", alpha=0.4)
plt.plot(np.repeat(data_bits, 1), 'go', label="True Bit (dots)", alpha=0.3)
plt.title(f"EMA Thresholds vs Log Power\nSNR = {SNR_dB} dB")
plt.xlabel("Bit Index")
plt.ylabel("Log Power / Threshold")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ==== PRINT BER RESULTS ====
print("📊 BER Results by α:")
for alpha, ber in ber_results.items():
    print(f"  α = {alpha:.2f} → BER = {ber:.2%}")
