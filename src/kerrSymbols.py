import sympy as sp
import numpy as np
import os
import pickle

CACHE_FILE = "gamma_expr_cache.pkl"

r, M, a = sp.symbols('r M a')  

if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, "rb") as f:
        Gamma_eq = pickle.load(f)
    print("Loaded cached Christoffel symbolic expressions.")
else:
    print("Computing Christoffel symbols symbolically (first run)...")
    t, theta, phi = sp.symbols('t theta phi')
    coords = [t, r, theta, phi]

    Sigma = r**2 + a**2*sp.cos(theta)**2
    Delta = r**2 - 2*M*r + a**2

    g = sp.zeros(4)
    g[0,0] = -(1 - 2*M*r / Sigma)
    g[1,1] = Sigma / Delta
    g[2,2] = Sigma
    g[3,3] = (r**2 + a**2 + 2*M*a**2*r*sp.sin(theta)**2 / Sigma) * sp.sin(theta)**2
    g[0,3] = -2*M*a*r*sp.sin(theta)**2 / Sigma
    g[3,0] = g[0,3]

    g_inv = g.inv()

    Gamma = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]

    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                expr = 0
                for sigma in range(4):
                    expr += g_inv[mu, sigma] * (
                        sp.diff(g[sigma, nu], coords[rho]) +
                        sp.diff(g[sigma, rho], coords[nu]) -
                        sp.diff(g[nu, rho], coords[sigma])
                    )
                Gamma[mu][nu][rho] = expr / 2 

    subs_eq = {theta: sp.pi/2, sp.cos(theta): 0, sp.sin(theta): 1}
    Gamma_eq = [[[Gamma[mu][nu][rho].subs(subs_eq) for rho in range(4)]
                 for nu in range(4)] for mu in range(4)]

    with open(CACHE_FILE, "wb") as f:
        pickle.dump(Gamma_eq, f)
    print("Cached Christoffel symbolic expressions.")

Gamma_funcs = [[[None for _ in range(4)] for _ in range(4)] for _ in range(4)]
for mu in range(4):
    for nu in range(4):
        for rho in range(4):
            expr = Gamma_eq[mu][nu][rho]
            if expr != 0:
                Gamma_funcs[mu][nu][rho] = sp.lambdify((r, M, a), expr, "numpy")

def kerr_christoffel_numeric(r_val, M_val, a_val):
    Γ = {}
    for mu in range(4):
        for nu in range(4):
            for rho in range(4):
                func = Gamma_funcs[mu][nu][rho]
                if func is not None:
                    val = func(r_val, M_val, a_val)
                    if abs(val) > 1e-12:
                        Γ[(mu, nu, rho)] = np.float64(val)
    return Γ