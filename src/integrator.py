import numpy as np
import math


def integrateGeodesic(geodesicFunc, state0, M, lambdaMax, stepInit=0.01, outputEvery=10):
    """
    Integrates a photon geodesic using adaptive RK4 with high precision.
    
    Parameters:
        geodesicFunc : function returning derivatives
        state0 : initial state [t, r, phi, tDot, rDot, phiDot]
        M : Schwarzschild mass
        lambdaMax : maximum affine parameter
        stepInit : initial step size
        outputEvery : record position every outputEvery steps

    Returns:
        positions : list of (x, y) tuples
        status : 'captured', 'escaped', or 'lambdaMaxReached'
    """
    
    state = np.array(state0, dtype=np.float64)
    lam = np.float64(0.0)
    step = np.float64(stepInit)
    M = np.float64(M)
    lambdaMax = np.float64(lambdaMax)
    outputEvery = int(outputEvery)

    positions = []
    stepCount = 0

    while lam < lambdaMax:
        t, r, phi, tDot, rDot, phiDot = state

        # --- Capture or escape conditions ---
        if r <= 2 * M:
            return positions, "captured"

        if r > 1e3 * M:
            return positions, "escaped"

        # Record every outputEvery steps
        if stepCount % outputEvery == 0:
            x = np.float64(r * math.cos(float(phi)))
            y = np.float64(r * math.sin(float(phi)))
            positions.append((float(x), float(y)))

        # --- Adaptive RK4 step ---
        k1 = np.array(geodesicFunc(state, M), dtype=np.float64)

        k2 = np.array(geodesicFunc(
            state + np.float64(0.5) * step * k1,
            M
        ), dtype=np.float64)

        k3 = np.array(geodesicFunc(
            state + np.float64(0.5) * step * k2,
            M
        ), dtype=np.float64)

        k4 = np.array(geodesicFunc(
            state + step * k3,
            M
        ), dtype=np.float64)

        newState = state + (step / np.float64(6.0)) * (k1 + 2 * k2 + 2 * k3 + k4)

        # --- Adaptive step sizing ---
        r_new = newState[1]
        
        if r_new < 4 * M:
            step = np.float64(np.maximum(0.001, float(step * 0.5)))
        elif r_new > 10 * M:
            step = np.float64(np.minimum(0.05, float(step * 1.2)))

        state = newState
        lam += step
        stepCount += 1

    return positions, "lambdaMaxReached"