import math

def integrateGeodesic(
    geodesicFunc,
    initialState,
    blackHoleMass,
    maxAffineParameter,
    initialStepSize=0.01,
    outputInterval=10,
    recordTrajectory=True,
    plotCallback=None
):
    state = initialState[:]
    affineParameter = 0.0
    stepSize = initialStepSize
    positions = []
    stepCounter = 0

    while affineParameter < maxAffineParameter:
        _, radius, angle, timeCoord, _, _ = state

        if radius <= 2*blackHoleMass:
            status = "captured"
            break
        if radius > 1e3*blackHoleMass:
            status = "escaped"
            break

        if recordTrajectory and stepCounter % outputInterval == 0:
            x = radius * math.cos(angle)
            y = radius * math.sin(angle)
            positions.append((x, y))
            if plotCallback:
                plotCallback(x, y)

        k1 = geodesicFunc(state, blackHoleMass)
        k2 = geodesicFunc([s + 0.5*stepSize*dk for s, dk in zip(state, k1)], blackHoleMass)
        k3 = geodesicFunc([s + 0.5*stepSize*dk for s, dk in zip(state, k2)], blackHoleMass)
        k4 = geodesicFunc([s + stepSize*dk for s, dk in zip(state, k3)], blackHoleMass)

        state = [
            s + (stepSize/6)*(k1_i + 2*k2_i + 2*k3_i + k4_i)
            for s, k1_i, k2_i, k3_i, k4_i in zip(state, k1, k2, k3, k4)
        ]

        if radius < 4*blackHoleMass:
            stepSize = max(0.001, stepSize*0.5)
        elif radius > 10*blackHoleMass:
            stepSize = min(0.05, stepSize*1.2)

        affineParameter += stepSize
        stepCounter += 1
    else:
        status = "escaped"

    finalTime = state[0]

    return positions, status, finalTime
