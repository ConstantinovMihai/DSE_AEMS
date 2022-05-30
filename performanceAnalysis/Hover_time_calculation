import math
import matplotlib.pyplot as plt

Ib = 14.4
Cb1 = 9500

Cbper_lst = []
Thover_lst = []
for Cb in range(0, Cb1, 1):
    Cmin = 0.2 * Cb1
    Thover = 0.06 * (Cb - Cmin) / Ib
    Thover_lst.append(Thover)
    Cbper = Cb*100/9500
    Cbper_lst.append(Cbper)

plt.plot(Cbper_lst, Thover_lst)
plt.xlabel('Battery Charge [%]')
plt.ylabel('Hover Time [min]')
plt.axis([0, 100, 0, 40])
plt.grid(True)
plt.show()
print("Hover time = ", Thover)
