# integrator.py

import math


def integrateGeodesic(geodesicFunc, state0, M, lambdaMax, stepInit=0.01, outputEvery=10):
    """
    Integrates a photon geodesic in Schwarzschild spacetime (equatorial plane).

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

    state = state0[:]
    lam = 0.0
    step = stepInit

    positions = []
    stepCount = 0

    while lam < lambdaMax:

        t, r, phi, tDot, rDot, phiDot = state

        # --- Capture or escape ---
        if r <= 2 * M:
            return positions, "captured"

        if r > 1e3 * M:
            return positions, "escaped"

        # Record every outputEvery steps
        if stepCount % outputEvery == 0:
            x = r * math.cos(phi)
            y = r * math.sin(phi)
            positions.append((x, y))

        # --- RK4 step ---
        k1 = geodesicFunc(state, M)

        k2 = geodesicFunc(
            [s + 0.5 * step * dk for s, dk in zip(state, k1)],
            M
        )

        k3 = geodesicFunc(
            [s + 0.5 * step * dk for s, dk in zip(state, k2)],
            M
        )

        k4 = geodesicFunc(
            [s + step * dk for s, dk in zip(state, k3)],
            M
        )

        newState = [
            s + (step / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
            for s, k1_i, k2_i, k3_i, k4_i in zip(state, k1, k2, k3, k4)
        ]

        # --- Adaptive step size ---
        if r < 4 * M:
            step = max(0.001, step * 0.5)

        elif r > 10 * M:
            step = min(0.05, step * 1.2)

        state = newState
        lam += step
        stepCount += 1

    return positions, "lambdaMaxReached"