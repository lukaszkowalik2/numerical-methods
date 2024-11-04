import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("times.csv")

plt.figure(figsize=(10, 6))
plt.plot(data["N"], data["Gauss"], label="Gauss")
plt.plot(data["N"], data["Thomas"], label="Thomas")
plt.plot(data["N"], data["LU"], label="LU")
plt.plot(data["N"], data["Eigen"], label="Eigen")

plt.xlabel("N")
plt.ylabel("Czas wykonania (s)")
plt.title("Zależność czasowa dla różnych algorytmów w funkcji N (co 10)")
plt.legend()
plt.grid()
plt.savefig("performance_chart.png")  
plt.show()
