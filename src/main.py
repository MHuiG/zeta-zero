import numpy as np
import sympy as sym
from scipy.special import gamma, comb

ACCURACY = 1e-10
ERROR = 1e-10
CHANGE_METHOD = 200

PI = np.pi
z = sym.Symbol('z')
Phi_0 = sym.cos(0.5 * sym.pi * (z ** 2) + 3 * sym.pi / 8) / sym.cos(sym.pi * z)
Phi_1 = (1 / (12 * sym.pi * sym.pi)) * sym.diff(Phi_0, z, 3)


def compute_zeta_Alternating_Series(t):

    mu = complex(0.5, -t)
    NUM = np.ceil((np.log2(1+2*t) + PI/2 * t * np.log2(np.e) - np.log2(ERROR) - np.log2(abs(1-2**mu))) / 3)

    def enum(k):
        i = k
        sub_sum = 0
        while i <= NUM:
            sub_sum += comb(NUM, i)
            i += 1

        return sub_sum

    zeta = complex(0, 0)
    s = complex(0.5, t)

    k = 1
    while k <= NUM:
        zeta += (-1)**(k-1) / k**s
        k += 1

    second = 0
    while k <= 2*NUM:
        second += (-1)**(k-1) * enum(k-NUM) / k**s
        k += 1

    second /= 2**NUM

    zeta += second
    zeta /= 1-2**(1-s)

    return zeta


f = open("./results/test.md", 'w')
f.write("|  No.   | Zero  |\n|  ----  | ----  |\n")

def write_zero(no,zero):
  f.write("|  "+str(no)+" | 1/2+"+str(zero)+"i |\n")


print("hello")
write_zero(1,11)
write_zero(2,22)
