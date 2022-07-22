import numpy as np
import matplotlib.pyplot as plt

Y, X = np.mgrid[-5:5:100j, -5:5:100j]
R = np.sqrt(X**2 + Y**2)
U = -1 - X**2 + Y
V = 1 + X - Y**2
U[R<1] = np.nan
V[R<1] = np.nan

plt.streamplot(X, Y, U, V, density=2.5, arrowstyle='-')
plt.axis("image")
plt.savefig("plots_to_sort/stream.png", dpi=300)
plt.show()
