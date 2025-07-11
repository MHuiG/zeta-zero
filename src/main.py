import numpy as np
import sympy as sym
from scipy.special import gamma, comb
import time

ACC=5
ACCURACY = 1e-5
ERROR = 1e-5
CHANGE_METHOD = 200 # 200 0 -1
ZeroNo = 0
ZeroHigh = 0
ZeroHighAdd = 1000 #1000
DELTA = 0.1

PI = np.pi
z = sym.Symbol('z')
Phi_0 = sym.cos(0.5 * sym.pi * (z ** 2) + 3 * sym.pi / 8) / sym.cos(sym.pi * z)
Phi_1 = (1 / (12 * sym.pi * sym.pi)) * sym.diff(Phi_0, z, 3)

p = sym.Symbol('p')
C_0 = sym.cos(2 * sym.pi * ((p ** 2) - p - (1 / 16))) / sym.cos(2 * sym.pi * p)
C_1 = - (1 / ((2 ** 5) * 3 * (sym.pi ** 2))) * sym.diff(C_0, p, 3)
C_2 = (1 / ((2 ** 11) * (3 ** 2) * (sym.pi ** 4))) * sym.diff(C_0, p, 6) + (1 / ((2 ** 6) * (sym.pi ** 2))) * sym.diff(C_0, p, 2)
C_3 = - (1 / ((2 ** 16) * (3 ** 4) * (sym.pi ** 6))) * sym.diff(C_0, p, 9) - (1 / ((2 ** 8) * 3 * 5 * (sym.pi ** 4))) * sym.diff(C_0, p, 5) - (1 / ((2 ** 6) * (sym.pi ** 2))) * sym.diff(C_0, p, 1)
C_4 = (1 / ((2 ** 23) * (3 ** 5) * (sym.pi ** 8))) * sym.diff(C_0, p, 12) + (11 / ((2 ** 17) * (3 ** 2) * 5 * (sym.pi ** 6))) * sym.diff(C_0, p, 8) + (19 / ((2 ** 13) * 3  * (sym.pi ** 4))) * sym.diff(C_0, p, 4) + (1 / ((2 ** 7) * (sym.pi ** 2))) * C_0

def compute_zeta_AS(t):
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

def compute_Zeta_AS(t):
    return (compute_zeta_AS(t) * np.exp(complex(0, compute_theta(t)))).real

def Chi(s):
    return 2**s * PI**(s-1) * np.sin(PI*s/2) * gamma(1-s)

def compute_theta(t):
    return (t / 2) * np.log(t / (2 * PI)) - (t / 2) - PI / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)

def compute_Phi(z0):
    return Phi_0.evalf(subs={z: z0}), Phi_1.evalf(subs={z: z0})

def compute_C(z0):
    return C_0.evalf(subs={p: z0}), C_1.evalf(subs={p: z0}), C_2.evalf(subs={p: z0}), C_3.evalf(subs={p: z0}), C_4.evalf(subs={p: z0})

def compute_Zeta_RS(t):
    tau = np.sqrt(t / (2 * PI))
    m = np.floor(tau)
    z0 = 2 * (tau - m) - 1
    Zeta = 0
    n = 1
    while n <= m:
        Zeta += 2 * (np.cos(compute_theta(t) - t * np.log(n))) / np.sqrt(n)
        n += 1
    phi = compute_Phi(z0)
    remain = phi[0] - phi[1] / tau
    remain *= (-1) ** (m + 1)
    remain /= np.sqrt(tau)
    Zeta += remain
    return Zeta

def compute_zeta_RS(t):
    mu = complex(0.5, -t)
    return np.exp(mu) * compute_Zeta_RS(t)

def compute_Zeta_RS_ACC(t):
    tau = np.sqrt(t / (2 * PI))
    m = np.floor(tau)
    z0 = tau - m
    Zeta = 0
    n = 1
    while n <= m:
        Zeta += 2 * (np.cos(compute_theta(t) - t * np.log(n))) / np.sqrt(n)
        n += 1
    phi = compute_C(z0)
    remain = phi[0] + (phi[1] / tau) + (phi[2] / (tau ** 2)) + (phi[3] / (tau ** 3)) + (phi[4] / (tau ** 4))
    remain *= (-1) ** (m - 1)
    remain /= np.sqrt(tau)
    Zeta += remain
    return Zeta

def compute_zeta_RS_ACC(t):
    mu = complex(0.5, -t)
    return np.exp(mu) * compute_Zeta_RS_ACC(t)

def zeros_numbers(T):
    estimate = compute_theta(T) / PI + 1
    if estimate - np.floor(estimate) < np.ceil(estimate) - estimate:
        return np.floor(estimate)
    else:
        return np.ceil(estimate)

def compute_Zeta(t):
    if CHANGE_METHOD == -1:
        return compute_Zeta_RS(t) #_ACC
    if (t < CHANGE_METHOD) and (t > 0):
        return compute_Zeta_AS(t)
    elif t >= CHANGE_METHOD:
        return compute_Zeta_RS(t)
    else:
        raise TypeError("Argument t should be a positive real number.")

def check_RH(T, delta):
    t1 = time.perf_counter()
    t = ZeroHigh
    count_zeros = ZeroNo
    while t < T:
        if np.sign(compute_Zeta(t)) * compute_Zeta(t + delta) < 0:
            count_zeros += 1
            print("Zero No.{}:\t({}, {})\n".format(count_zeros, np.round(t, ACC), np.round(t+delta, ACC)))
            compute(count_zeros,np.round(t, ACC),np.round(t+delta, ACC))
        t += delta
    print("Find {} zeros with 0<t<{}.".format(count_zeros, T))
    print("Expecting {} zeros.".format(zeros_numbers(T)))
    print("\n")
    print("Expecting {} zeros.\n".format(zeros_numbers(T)))
    print("Find {} zeros.\n".format(count_zeros))
    t2 = time.perf_counter()
    print("Total time cost: {} seconds.".format(t2 - t1))
    print("Average time cost: {} seconds per zero.".format(round((t2 - t1) / count_zeros, 7)))
    f_index.write(str(ZeroHigh+ZeroHighAdd)+"\n")
    f_index.write(str(count_zeros)+"\n")
    if count_zeros == zeros_numbers(T):
        return True
    else:
        return False

def compute_zero(low, high, method, accuracy=ACCURACY):
    if np.sign(method(low)) * np.sign(method(high)) > 0:
        raise ValueError("Suspect there is no zero between {} and {}! Please check the interval.".format(low, high))
    else:
        err = high - low
        mid = (high + low) / 2
        low_value = method(low)
        mid_value = method(mid)
        high_value = method(high)
        while err > accuracy:
            if mid_value == 0:
                return mid
            elif np.sign(mid_value) * np.sign(low_value) < 0:
                high = mid
            elif np.sign(mid_value) * np.sign(high_value) < 0:
                low = mid
            err = high - low
            mid = (high + low) / 2
            low_value = method(low)
            mid_value = method(mid)
            high_value = method(high)
        digit = -int(np.log10(accuracy))
        return '%.{}f'.format(digit) % mid

def compute(num,low,high):
    if num < CHANGE_METHOD:
        method = compute_Zeta_AS
    else:
        method = compute_Zeta_RS
    if CHANGE_METHOD == -1:
        method = compute_Zeta_RS_ACC
    zero = compute_zero(low, high, method)
    print("Zero No {}:\t".format(num))
    print(str(zero))
    print("\n")
    write_zero(num,zero)

if __name__ == "__main__":
    with open('./results/index.md', 'r') as file: 
        content = file.read() 
        print(content)
        a = content.split()
        ZeroHigh = int(a[0])
        ZeroNo = int(a[1])
        file.close()

    f = open("./results/"+str(ZeroHigh)+"-"+str(ZeroHigh+ZeroHighAdd)+".md", 'w')
    f.write("|  No.   | Zero  |\n|  ----  | ----  |\n")

    f_index = open("./results/index.md", 'w')

    def write_zero(no,zero):
      f.write("|  "+str(no)+" | 1/2+"+str(zero)+"i |\n")

    
    res=check_RH(ZeroHigh+ZeroHighAdd, DELTA)
    print(res)
"""
    if not res:
        fe = open("./results/errors/error.md", 'w')
        fe.write("RH is not valid.\nHight:"+str(ZeroHigh+ZeroHighAdd))
"""
