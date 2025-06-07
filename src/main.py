import numpy as np
import sympy as sym

PI = np.pi
z = sym.Symbol('z')
Phi_0 = sym.cos(0.5 * sym.pi * (z ** 2) + 3 * sym.pi / 8) / sym.cos(sym.pi * z)
Phi_1 = (1 / (12 * sym.pi * sym.pi)) * sym.diff(Phi_0, z, 3)





f = open("./results/test.md", 'w')
f.write("|  No.   | Zero  |\n|  ----  | ----  |\n")

def write_zero(no,zero):
  f.write("|  "+str(no)+" | 1/2+"+str(zero)+"i |\n")


print("hello")
write_zero(1,11)
write_zero(2,22)
