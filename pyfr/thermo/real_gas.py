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


def rho_from_pT(p, T, R=R, a=A, b=B, tol=1e-10, maxiter=20):
    """Return the density from pressure and temperature via VdW EOS."""
    # Initial guess from the ideal-gas relation
    rho = p / (R * T)

    for _ in range(maxiter):
        f = pressure(rho, T, R, a, b) - p
        df = (R * T) / ((1 - b * rho) ** 2) - 2 * a * rho
        rho_new = rho - f / df

        if abs(rho_new - rho) < tol:
            return rho_new
        rho = rho_new

    return rho


def sound_speed_sq(rho, T, R=R, a=A, b=B, cv=CV):
    """Square of the adiabatic sound speed for a Van der Waals gas."""
    ombr = 1.0 - b * rho
    return (R * T / (ombr * ombr)) * (1.0 + R / cv) - 2.0 * a * rho
