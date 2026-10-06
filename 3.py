import math

def f(x, t):
    return math.cos(t / (1.0 + x * x) + 0.001 * x)

def trapezoid(t, n, a, b):
    h = (b - a) / n
    s = 0.5 * (f(a, t) + f(b, t))
    for i in range(1, n):
        s += f(a + i * h, t)
    return h * s

def runge_rule(t, eps, a, b):
    n = 1
    I_old = trapezoid(t, n, a, b)
    iter_count = 1
    while True:
        n *= 2
        I_new = trapezoid(t, n, a, b)
        iter_count += 1
        if abs(I_new - I_old) < eps or n > 1048576:
            return I_new, iter_count
        I_old = I_new

def gauss3(t, a, b):
    x3 = [-0.7745966692414834, 0.0, 0.7745966692414834]
    w3 = [5.0/9.0, 8.0/9.0, 5.0/9.0]
    tr = (b - a) / 2
    sh = (a + b) / 2
    sum_val = 0
    for i in range(3):
        sum_val += w3[i] * f(tr * x3[i] + sh, t)
    return tr * sum_val

def gauss4(t, a, b):
    x4 = [-0.8611363115940526, -0.3399810435848563,
          0.3399810435848563,  0.8611363115940526]
    w4 = [0.3478548451374538, 0.6521451548625461,
          0.6521451548625461, 0.3478548451374538]
    tr = (b - a) / 2
    sh = (a + b) / 2
    sum_val = 0
    for i in range(4):
        sum_val += w4[i] * f(tr * x4[i] + sh, t)
    return tr * sum_val

def main():
    a = 0.0
    b = 2.0
    c = 0.0
    d = 1.0
    m = 20
    eps = 0.001
    r = (d - c) / m

    print(f"{'j':<3}{'t':<12}{'Runge (trapezoid)':<20}{'Gauss3':<16}{'Gauss4':<16}{'iter':<6}")
    print("-" * 75)

    for j in range(m + 1):
        t = c + j * r
        I_runge, iter_count = runge_rule(t, eps, a, b)
        I_g3 = gauss3(t, a, b)
        I_g4 = gauss4(t, a, b)
        print(f"{j:<3}{t:<12.8f}{I_runge:<20.8f}{I_g3:<16.8f}{I_g4:<16.8f}{iter_count:<6}")

if __name__ == "__main__":
    main()
