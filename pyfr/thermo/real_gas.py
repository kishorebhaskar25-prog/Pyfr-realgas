R = 287.0
A = 0.9
B = 0.1
CV = 718.0

def pressure(rho, T, R=R, a=A, b=B):
    v = 1.0 / rho
    return (R * T) / (v - b) - (a / (v ** 2))

def energy(rho, T, cv=CV, a=A):
    return cv * T - a * rho

def T_from_rho_e(rho, e, cv=CV, a=A):
    return (e + a * rho) / cv
