import numpy as np

def schwarzschildMetric(r, M):
    """
    Return the Schwarzschild metric tensor in equatorial plane (θ = π/2).
    High precision version using numpy float64.
    """
    r = np.float64(r)
    M = np.float64(M)
    
    g = np.zeros((4, 4), dtype=np.float64)
    
    g[0, 0] = -(1 - 2 * M / r)
    g[1, 1] = 1 / (1 - 2 * M / r)
    g[2, 2] = r**2
    g[3, 3] = r**2

    return g


def christoffelSymbols(r, M):
    """
    Return non-zero Christoffel symbols for equatorial motion.
    High precision version using numpy float64.
    """
    r = np.float64(r)
    M = np.float64(M)
    
    Γ = {}

    Γ[(0, 1, 0)] = M / (r * (r - 2 * M))
    Γ[(1, 0, 0)] = M * (r - 2 * M) / r**3
    Γ[(1, 1, 1)] = -M / (r * (r - 2 * M))
    Γ[(1, 3, 3)] = -(r - 2 * M)
    Γ[(3, 1, 3)] = 1 / r
    Γ[(3, 3, 1)] = 1 / r

    return Γ


def geodesicEquations(state, M):
    """
    state: [t, r, phi, tDot, rDot, phiDot]
    
    returns: [tDot, rDot, phiDot, tDDot, rDDot, phiDDot]
    
    High precision version using numpy float64.
    """
    state = np.array(state, dtype=np.float64)
    M = np.float64(M)
    
    t, r, phi, tDot, rDot, phiDot = state

    Γ = christoffelSymbols(r, M)

    tDDot = -2 * Γ.get((0, 1, 0), np.float64(0)) * tDot * rDot

    rDDot = (
        -Γ.get((1, 0, 0), np.float64(0)) * tDot**2
        - Γ.get((1, 1, 1), np.float64(0)) * rDot**2
        - Γ.get((1, 3, 3), np.float64(0)) * phiDot**2
    )

    phiDDot = -2 * Γ.get((3, 1, 3), np.float64(0)) * rDot * phiDot

    return np.array([tDot, rDot, phiDot, tDDot, rDDot, phiDDot], dtype=np.float64)