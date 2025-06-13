import usb.core
import usb.util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# === RTL-SDR Parameters ===
VENDOR_ID = 0x0bda
PRODUCT_ID = 0x2832
EP_IN = 0x81
IFACE = 0
FRAME_SIZE = 512 * 128  # = 65536 bytes

# === Frequency ที่ต้องการรับ ===
target_freq_hz = 100_000_000  # 🧠 ตั้งตรงนี้ เช่น 100e6 หรือ 433920000

# === หาอุปกรณ์ ===
dev = usb.core.find(idVendor=VENDOR_ID, idProduct=PRODUCT_ID)
if dev is None:
    raise ValueError("❌ ไม่พบ RTL2832U")

try:
    if dev.is_kernel_driver_active(IFACE):
        dev.detach_kernel_driver(IFACE)
except (NotImplementedError, usb.core.USBError):
    pass

usb.util.claim_interface(dev, IFACE)

# === ส่ง Vendor Command ===
def ctrl(index, value, data=None):
    return dev.ctrl_transfer(0x40, index, value, 0, data or [])

# ✅ เริ่มต้นอุปกรณ์
# ctrl(0x01, 0x00)  # อย่าส่ง reset (เสี่ยง I/O Error)

# 1. ตั้ง sample rate = 2.048 MSPS
ctrl(0x1e, 0x0001, [0x58, 0x1b])  # = 0x1B58 → 7000 = 2.048MSPS

# 2. ตั้งความถี่ = target_freq_hz (LSB first 4 byte)
freq_bytes = target_freq_hz.to_bytes(4, byteorder='little')
ctrl(0x19, 0x00, list(freq_bytes))

# 3. เปิด stream
ctrl(0x15, 0x00)

# === วาดกราฟ Spectrum
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=1)
ax.set_ylim(-80, 0)
ax.set_xlim(-1, 1)
ax.set_title(f"Live Spectrum at {target_freq_hz/1e6:.2f} MHz")
ax.set_xlabel("Normalized Frequency")
ax.set_ylabel("Power (dB)")

def update(_):
    try:
        raw = dev.read(EP_IN, FRAME_SIZE, timeout=1000)
        data = np.frombuffer(raw, dtype=np.uint8)
        I = data[::2] - 127.5
        Q = data[1::2] - 127.5
        iq = I + 1j * Q
        spectrum = np.fft.fftshift(np.fft.fft(iq))
        power = 10 * np.log10(np.abs(spectrum)**2 + 1e-12)
        x = np.linspace(-1, 1, len(power))
        line.set_data(x, power)
        return line,
    except usb.core.USBError as e:
        print(f"⚠️ USBError: {e}")
        return line,

ani = animation.FuncAnimation(fig, update, interval=500, blit=True)
plt.show()

# === ปิดการใช้งาน
usb.util.release_interface(dev, IFACE)
dev.reset()
