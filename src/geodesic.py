def schwarzschildMetric(radius, blackHoleMass):
    g = [[0]*4 for _ in range(4)]
    g[0][0] = -(1 - 2*blackHoleMass/radius)
    g[1][1] = 1/(1 - 2*blackHoleMass/radius)
    g[2][2] = radius**2
    g[3][3] = radius**2
    return g

def christoffelSymbols(radius, blackHoleMass):
    christoffel = {}
    christoffel[(0,1,0)] = blackHoleMass / (radius*(radius - 2*blackHoleMass))
    christoffel[(1,0,0)] = blackHoleMass*(radius - 2*blackHoleMass)/radius**3
    christoffel[(1,1,1)] = -blackHoleMass/(radius*(radius - 2*blackHoleMass))
    christoffel[(1,3,3)] = -(radius - 2*blackHoleMass)
    christoffel[(3,1,3)] = 1/radius
    return christoffel

def geodesicEquations(state, blackHoleMass):
    timeCoord, radius, angle, timeVelocity, radialVelocity, angularVelocity = state
    christoffel = christoffelSymbols(radius, blackHoleMass)

    timeAcceleration = -2 * christoffel.get((0,1,0),0) * timeVelocity * radialVelocity
    radialAcceleration = -christoffel.get((1,0,0),0)*timeVelocity**2 - christoffel.get((1,1,1),0)*radialVelocity**2 - christoffel.get((1,3,3),0)*angularVelocity**2
    angularAcceleration = -2 * christoffel.get((3,1,3),0) * radialVelocity * angularVelocity

    return [timeVelocity, radialVelocity, angularVelocity, timeAcceleration, radialAcceleration, angularAcceleration]
