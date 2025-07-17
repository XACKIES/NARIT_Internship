import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal

# ==== PARAMETERS ====
bit_rate = 1_000_000
fs = 42_000_000
samples_per_bit = fs // bit_rate
num_bits = 50000
epsilon = 1e-12
np.random.seed(42)

SNR_dB_range = np.arange(0, 11, 1)
alpha_values = [0.1,0.3, 0.5, 0.7, 0.9]
ber_results = {}

# ==== LOOP: EMA Threshold ====
for alpha in alpha_values:
    ber_list = []
    for SNR_dB in SNR_dB_range:
        data_bits = np.random.randint(0, 2, num_bits)
        ook_signal = np.repeat(data_bits, samples_per_bit).astype(float)

        signal_power = np.mean(ook_signal ** 2)
        SNR_linear = 10 ** (SNR_dB / 10)
        noise_power = signal_power / SNR_linear
        noise = normal(0, np.sqrt(noise_power), len(ook_signal))
        received_signal = ook_signal + noise

        log_power = []
        for i in range(num_bits):
            segment = received_signal[i * samples_per_bit : (i + 1) * samples_per_bit]
            energy = np.sum(segment ** 2)
            log_power.append(np.log10(energy + epsilon))
        log_power = np.array(log_power)

        threshold = np.zeros(num_bits)
        threshold[0] = log_power[0]
        for i in range(1, num_bits):
            threshold[i] = alpha * log_power[i] + (1 - alpha) * threshold[i - 1]

        detected_bits = (log_power > threshold).astype(int)
        ber = np.sum(detected_bits != data_bits) / num_bits
        ber_list.append(ber)
    ber_results[f"EMA α={alpha}"] = ber_list
# === Fixed Point Threshold (จากค่าเฉลี่ย log power ที่ SNR = 5 dB) ===
SNR_ref_dB = 10
data_bits_ref = np.random.randint(0, 2, num_bits)
ook_signal_ref = np.repeat(data_bits_ref, samples_per_bit).astype(float)
signal_power_ref = np.mean(ook_signal_ref ** 2)
SNR_ref_linear = 10 ** (SNR_ref_dB / 10)
noise_power_ref = signal_power_ref / SNR_ref_linear
noise_ref = normal(0, np.sqrt(noise_power_ref), len(ook_signal_ref))
received_signal_ref = ook_signal_ref + noise_ref

log_power_ref = []
for i in range(num_bits):
    segment = received_signal_ref[i * samples_per_bit : (i + 1) * samples_per_bit]
    energy = np.sum(segment ** 2)
    log_power_ref.append(np.log10(energy + epsilon))
fixed_point_threshold = np.mean(log_power_ref)

# === คำนวณ BER จาก Fixed Point Threshold ที่คงที่ตลอด ===
ber_fixed_point = []
for SNR_dB in SNR_dB_range:
    data_bits = np.random.randint(0, 2, num_bits)
    ook_signal = np.repeat(data_bits, samples_per_bit).astype(float)
    signal_power = np.mean(ook_signal ** 2)
    SNR_linear = 10 ** (SNR_dB / 10)
    noise_power = signal_power / SNR_linear
    noise = normal(0, np.sqrt(noise_power), len(ook_signal))
    received_signal = ook_signal + noise

    log_power = []
    for i in range(num_bits):
        segment = received_signal[i * samples_per_bit : (i + 1) * samples_per_bit]
        energy = np.sum(segment ** 2)
        log_power.append(np.log10(energy + epsilon))
    log_power = np.array(log_power)

    detected_bits = (log_power > fixed_point_threshold).astype(int)
    ber = np.sum(detected_bits != data_bits) / num_bits
    ber_fixed_point.append(ber)

ber_results["Fixed Point (from SNR=5dB)"] = ber_fixed_point


# ==== PLOT ====
plt.figure(figsize=(10, 6))
for label, ber_list in ber_results.items():
    plt.plot(SNR_dB_range, ber_list, marker='o', label=label)

plt.title("BER vs SNR for OOK Detection\n(EMA vs Fixed Threshold Recomputed per SNR)")
plt.xlabel("SNR (dB)")
plt.ylabel("Bit Error Rate (BER)")
plt.grid(True)
plt.legend()
plt.yscale("log")
plt.tight_layout()
plt.show()
