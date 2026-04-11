import numpy as np

def schwarzschildChristoffel(r, M):

    r = np.float64(r)
    M = np.float64(M)
    
    Γ = {}
    Γ[(0, 1, 0)] = M / (r * (r - 2 * M))
    Γ[(1, 0, 0)] = M * (r - 2 * M) / r**3
    Γ[(1, 1, 1)] = -M / (r * (r - 2 * M))
    Γ[(1, 3, 3)] = -(r - 2 * M)
    Γ[(3, 1, 3)] = 1 / r
    return Γ


from kerrSymbols import kerr_christoffel_numeric

def kerrChristoffel(r, M, a):
    return kerr_christoffel_numeric(r, M, a)


def geodesicEquations(state, M, a=0.0):

    state = np.array(state, dtype=np.float64)
    M = np.float64(M)
    a = np.float64(a)
    
    t, r, phi, tDot, rDot, phiDot = state
    
    if abs(a) < 1e-14:
        Γ = schwarzschildChristoffel(r, M)
    else:
        Γ = kerrChristoffel(r, M, a)
    
 
    tDDot = (
        -2 * Γ.get((0, 1, 0), np.float64(0)) * tDot * rDot
        - 2 * Γ.get((0, 1, 3), np.float64(0)) * rDot * phiDot
    )
    

    rDDot = (
        -Γ.get((1, 0, 0), np.float64(0)) * tDot**2
        - Γ.get((1, 1, 1), np.float64(0)) * rDot**2
        - Γ.get((1, 3, 3), np.float64(0)) * phiDot**2
    )
    

    phiDDot = (
        -2 * Γ.get((3, 1, 0), np.float64(0)) * tDot * rDot
        - 2 * Γ.get((3, 1, 3), np.float64(0)) * rDot * phiDot
    )
    
    return np.array([tDot, rDot, phiDot, tDDot, rDDot, phiDDot], dtype=np.float64)