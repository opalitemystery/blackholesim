import math

from schwarzchild import Gravity, Lightspeed, Mass, schwarzchildFactor, timeRate, angleRate

def radialDerivative(E, L, radius):
    return -(math.sqrt(E**2 - schwarzchildFactor(radius) * L**2 / radius**2))

def geodesicRHS(path_parameter, state, parameters):
    t, r, angle = state
    E = parameters["E"]
    L = parameters["L"]
    
    dt = timeRate(E, r)
    dr = radialDerivative(E, L, r)
    dangle = angleRate(L, r)
    
    return [dt, dr, dangle]