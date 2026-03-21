from geodesic import geodesicEquations
from integrator import integrateGeodesic
import numpy as np


def computeNullCondition(M, r0, rDot0, phiDot0):
    """
    Compute tDot0 from null geodesic condition.
    For a photon: g_μν dx^μ dx^ν = 0
    Uses high precision numpy float64
    """
    M = np.float64(M)
    r0 = np.float64(r0)
    rDot0 = np.float64(rDot0)
    phiDot0 = np.float64(phiDot0)
    
    # Photon must start outside event horizon (r0 > 2M)
    if r0 <= 2 * M:
        raise ValueError(f"Photon starting radius (r0={r0}) must be outside event horizon (2M={2*M})")
    
    factor = 1 - 2 * M / r0
    tDot0_squared = (rDot0**2 / factor + r0**2 * phiDot0**2) / factor
    if tDot0_squared < 0:
        raise ValueError("Invalid initial conditions: negative tDot0_squared")
    return np.sqrt(tDot0_squared)


def runSingleSimulation(M, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery):
    """
    Run a single photon geodesic simulation with high precision.
    
    Returns:
        trajectory: list of {"x": float, "y": float} dicts
        status: "captured", "escaped", or "lambdaMaxReached"
    """
    # Convert all inputs to float64 for high precision
    M = np.float64(M)
    r0 = np.float64(r0)
    phi0 = np.float64(phi0)
    rDot0 = np.float64(rDot0)
    phiDot0 = np.float64(phiDot0)
    lambdaMax = np.float64(lambdaMax)
    stepInit = np.float64(stepInit)
    outputEvery = int(outputEvery)
    
    tDot0 = computeNullCondition(M, r0, rDot0, phiDot0)
    state0 = np.array([np.float64(0.0), r0, phi0, tDot0, rDot0, phiDot0], dtype=np.float64)
    
    positions, status = integrateGeodesic(
        geodesicEquations,
        state0,
        M,
        lambdaMax,
        stepInit,
        outputEvery
    )
    
    trajectory = [{"x": float(x), "y": float(y)} for x, y in positions]
    return trajectory, status


# --- Optional: CLI interface for testing ---
if __name__ == "__main__":
    inputStr = input(
        "Enter the following values separated by spaces:\n"
        "M r0 phi0 rDot0 phiDot0 lambdaMax stepInit outputEvery\n"
        "Example: 1.0 5.0 0.0 0.0 0.192 100 0.01 10\n"
    )
    
    inputs = inputStr.strip().split()
    
    if len(inputs) != 8:
        raise ValueError("Expected 8 inputs")
    
    M = np.float64(inputs[0])
    r0 = np.float64(inputs[1])
    phi0 = np.float64(inputs[2])
    rDot0 = np.float64(inputs[3])
    phiDot0 = np.float64(inputs[4])
    lambdaMax = np.float64(inputs[5])
    stepInit = np.float64(inputs[6])
    outputEvery = int(inputs[7])
    
    trajectory, status = runSingleSimulation(
        M, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery
    )
    
    print("\nSimulation status:", status)
    print("Trajectory points:", len(trajectory))
    print("First few points:", trajectory[:3])