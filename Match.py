import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve

# 1. สร้างสัญญาณต้นแบบ s(t)
T = 50  # ความยาวสัญญาณ
s = np.ones(T)  # ใช้ rectangular pulse เป็นตัวอย่าง

# 2. สร้างสัญญาณที่มี noise
N = 200  # ความยาวของสัญญาณรับ
r = np.random.normal(0, 1, N)  # white Gaussian noise

# ฝังสัญญาณ s(t) ไว้ใน noise ที่ตำแหน่งที่รู้ (เช่น t = 80)
r[80:80+T] += s

# 3. สร้าง matched filter (time-reversed ของ s)
h = s[::-1]  # reversed s(t)

# 4. ทำ convolution (matched filtering)
y = convolve(r, h, mode='same')  # matched filter output

# 5. Plot ผลลัพธ์
plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(s)
plt.title("Template Signal s(t)")

plt.subplot(3, 1, 2)
plt.plot(r)
plt.title("Received Signal r(t) = s(t) + noise")

plt.subplot(3, 1, 3)
plt.plot(y)
plt.title("Matched Filter Output y(t)")
plt.xlabel("Time")

plt.tight_layout()
plt.show()
