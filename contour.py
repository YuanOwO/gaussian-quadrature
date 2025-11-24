import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(8, 6), tight_layout=True)
ax = fig.add_subplot(111)


def f(x, y):
    return 3 * np.sin(8 * np.pi * x) * np.cos(8 * np.pi * y) + x + y + 1


num = 4000
x = np.linspace(2, 6, num)
y = np.linspace(2, 6, num)
X, Y = np.meshgrid(x, y)
Z = f(X, Y)
contour = ax.contourf(X, Y, Z, levels=50, cmap="viridis")
fig.colorbar(contour)

ax.set_title("$f(x, y) = 3 \\sin(8 \\pi x) \\cos(8 \\pi y) + x + y + 1$")

fig.savefig("image/contour_plot.png", dpi=300)
