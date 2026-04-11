from geodesic import geodesicEquations
from integrator import integrateGeodesic
import numpy as np


def computeNullCondition(M, a, r0, rDot0, phiDot0):

    M = np.float64(M)
    a = np.float64(a)
    r0 = np.float64(r0)
    rDot0 = np.float64(rDot0)
    phiDot0 = np.float64(phiDot0)
    
    r_plus = M + np.sqrt(M*M - a*a)
    if r0 <= r_plus:
        raise ValueError(f"Photon starting radius (r0={r0}) must be outside event horizon (r+={r_plus})")
    

    r2 = r0 * r0
    a2 = a * a
    sigma = r2
    delta = r2 - 2*M*r0 + a2
    
    g_tt = -(1 - 2*M*r0 / sigma)
    g_rr = sigma / delta
    g_pp = r2 + a2 + 2*M*a2*r0 / sigma
    g_tp = -2*M*a*r0 / sigma
    
  
    
    A = g_tt
    B = 2 * g_tp * phiDot0
    C = g_rr * rDot0**2 + g_pp * phiDot0**2
    
    discriminant = B*B - 4*A*C
    if discriminant < 0:
        raise ValueError("Invalid initial conditions: negative discriminant")
    
    tDot0 = (-B + np.sqrt(discriminant)) / (2*A)
    
    return tDot0


def runSingleSimulation(M, a, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery):

    M = np.float64(M)
    a = np.float64(a)
    r0 = np.float64(r0)
    phi0 = np.float64(phi0)
    rDot0 = np.float64(rDot0)
    phiDot0 = np.float64(phiDot0)
    lambdaMax = np.float64(lambdaMax)
    stepInit = np.float64(stepInit)
    outputEvery = int(outputEvery)
    
    tDot0 = computeNullCondition(M, a, r0, rDot0, phiDot0)
    state0 = np.array([np.float64(0.0), r0, phi0, tDot0, rDot0, phiDot0], dtype=np.float64)
    
    def kerr_geodesic(state, M_val, a_val):
        return geodesicEquations(state, M_val, a_val)
    
    positions, status = integrateGeodesic(
        kerr_geodesic,
        state0,
        M,
        a, 
        lambdaMax,
        stepInit,
        outputEvery
    )
    
    trajectory = [{"x": float(x), "y": float(y)} for x, y in positions]
    return trajectory, status


if __name__ == "__main__":
    inputStr = input(
        "Enter the following values separated by spaces:\n"
        "M a r0 phi0 rDot0 phiDot0 lambdaMax stepInit outputEvery\n"
        "Example: 1.0 0.5 5.0 0.0 0.0 0.192 100 0.01 10\n"
    )
    
    inputs = inputStr.strip().split()
    
    if len(inputs) != 9:
        raise ValueError("Expected 9 inputs")
    
    M = np.float64(inputs[0])
    a = np.float64(inputs[1])
    r0 = np.float64(inputs[2])
    phi0 = np.float64(inputs[3])
    rDot0 = np.float64(inputs[4])
    phiDot0 = np.float64(inputs[5])
    lambdaMax = np.float64(inputs[6])
    stepInit = np.float64(inputs[7])
    outputEvery = int(inputs[8])
    
    trajectory, status = runSingleSimulation(
        M, a, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery
    )
    
    print("\nSimulation status:", status)
    print("Trajectory points:", len(trajectory))
    print("First few points:", trajectory[:3])