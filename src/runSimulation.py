from geodesic import geodesicEquations
from integrator import integrateGeodesic
import math


def computeNullCondition(M, r0, rDot0, phiDot0):
    """
    Compute tDot0 from null geodesic condition.
    For a photon: g_μν dx^μ dx^ν = 0
    """
    factor = 1 - 2 * M / r0
    tDot0_squared = (rDot0**2 / factor + r0**2 * phiDot0**2) / factor
    if tDot0_squared < 0:
        raise ValueError("Invalid initial conditions")
    return math.sqrt(tDot0_squared)


def runSingleSimulation(M, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery):
    """
    Run a single photon geodesic simulation.
    
    Returns:
        trajectory: list of {"x": float, "y": float} dicts
        status: "captured", "escaped", or "lambdaMaxReached"
    """
    tDot0 = computeNullCondition(M, r0, rDot0, phiDot0)
    state0 = [0.0, r0, phi0, tDot0, rDot0, phiDot0]
    
    positions, status = integrateGeodesic(
        geodesicEquations,
        state0,
        M,
        lambdaMax,
        stepInit,
        outputEvery
    )
    
    trajectory = [{"x": x, "y": y} for x, y in positions]
    return trajectory, status


def runBatchSimulation(M, r0, phi0, rDot0, phiDot0Array, lambdaMax, stepInit, outputEvery):
    """
    Run multiple photon geodesic simulations with different angular velocities.
    
    Returns:
        list of dicts with keys: phiDot0, status, trajectory
    """
    results = []
    
    for phiDot0 in phiDot0Array:
        try:
            trajectory, status = runSingleSimulation(
                M, r0, phi0, rDot0, float(phiDot0), lambdaMax, stepInit, outputEvery
            )
            results.append({
                "phiDot0": float(phiDot0),
                "status": status,
                "trajectory": trajectory
            })
        except Exception as e:
            results.append({
                "phiDot0": float(phiDot0),
                "error": str(e)
            })
    
    return results


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
    
    M = float(inputs[0])
    r0 = float(inputs[1])
    phi0 = float(inputs[2])
    rDot0 = float(inputs[3])
    phiDot0 = float(inputs[4])
    lambdaMax = float(inputs[5])
    stepInit = float(inputs[6])
    outputEvery = int(inputs[7])
    
    trajectory, status = runSingleSimulation(
        M, r0, phi0, rDot0, phiDot0, lambdaMax, stepInit, outputEvery
    )
    
    print("\nSimulation status:", status)
    print("Trajectory points:", len(trajectory))
    print("First few points:", trajectory[:3])