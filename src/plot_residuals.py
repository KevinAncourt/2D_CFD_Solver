import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../output/residual_history.dat", comments="#")

iterations = data[:, 0]
residuals  = data[:, 1]

plt.figure()
plt.loglog(iterations, residuals)
plt.xlabel("Iteration")
plt.ylabel("L2 residual (rho)")
plt.title("Residual convergence history")
plt.grid(True)
plt.tight_layout()
plt.show()