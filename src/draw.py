import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from main import compute_zeta_AS

def draw_zeta_critical_line(begin, end, step=0.1):
    fig = plt.figure(figsize=(10, 10))
    plt.xlabel(r'$\mathfrak{R}\zeta(\frac{1}{2}+it)$')
    plt.ylabel(r'$\mathfrak{I}\zeta(\frac{1}{2}+it)$')
    plt.xlim(-2, 4)
    plt.ylim(-3, 3)
    plt.grid(linestyle=':')
    plt.title('Dynamical '+r'$\zeta(s)$'+' on the critical line')
    x, y = [], []

    def update(n):
        x.append(compute_zeta_AS(n).real)
        y.append(compute_zeta_AS(n).imag)
        plt.plot(x, y, "r-", linewidth=0.5, label=r'$\zeta(\frac{1}{2}+it)$')

    ani = FuncAnimation(fig, update, frames=np.arange(begin, end, step), interval=50, blit=False, repeat=False)
    ani.save("./results/dynamical_zeta.gif", writer='pillow')
    #plt.show()

    return


if __name__ == "__main__":
    draw_zeta_critical_line(1, 100)

    
