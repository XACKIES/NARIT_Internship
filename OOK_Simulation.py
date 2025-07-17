import numpy as np
import matplotlib.pyplot as plt

# === PARAMETERS ===
Fs = 500e6              # Sampling rate: 500 MHz
Fc = 1090e6             # Carrier frequency: 1090 MHz
LO_freq = 1040e6        # LO frequency: 1 GHz
IF = Fc - LO_freq       # Intermediate frequency: 90 MHz
bitrate = 1e6           # Bitrate: 1 Mbps → 1 bit = 1 µs
bit_pattern = [1, 0,1,0,1,0] * 1  
num_bits = len(bit_pattern)
duration = num_bits / bitrate  

# === TIME VECTOR ===
t = np.linspace(0, duration, int(Fs * duration), endpoint=False)

# === OOK MODULATION ===
samples_per_bit = int(Fs / bitrate)
symbols = np.repeat(bit_pattern, samples_per_bit)
ook_signal = symbols[:len(t)]

carrier = np.cos(2 * np.pi * Fc * t)
modulated = ook_signal * carrier

# === DOWNCONVERSION ===
lo = np.cos(2 * np.pi * LO_freq * t)
downconverted = modulated * lo

# === LOW-PASS FILTER ===
def moving_average(x, N):
    if N <= 0:
        raise ValueError("Filter window size N must be > 0")
    return np.convolve(x, np.ones(N)/N, mode='same')

N = max(1, int(Fs / (2 * 200e6)))  
filtered = moving_average(downconverted, N)

# === PLOTS ===
plt.figure(figsize=(14, 9))

plt.subplot(3, 1, 1)
plt.plot(t * 1e6, modulated)
plt.title(f'OOK Modulated Signal @ 1090 MHz (Pattern "1010" x 25)')
plt.xlabel('Time (µs)')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 2)
plt.plot(t * 1e6, downconverted)
plt.title(f'Downconverted Signal (IF = {IF/1e6:.1f} MHz)')
plt.xlabel('Time (µs)')
plt.ylabel('Amplitude')

plt.subplot(3, 1, 3)
plt.plot(t * 1e6, filtered)
plt.title('Filtered Signal (Moving Average LPF)')
plt.xlabel('Time (µs)')
plt.ylabel('Amplitude')

plt.tight_layout()
plt.show()
