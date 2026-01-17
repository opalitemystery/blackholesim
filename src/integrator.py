import math

def rk4Step(f, pathParameter, state, stepSize, parameters):
  
    k1 = f(pathParameter, state, parameters)
    k2 = f(pathParameter + stepSize / 2, [s + stepSize / 2 * k for s, k in zip(state, k1)], parameters)
    k3 = f(pathParameter + stepSize / 2, [s + stepSize / 2 * k for s, k in zip(state, k2)], parameters)
    k4 = f(pathParameter + stepSize, [s + stepSize * k for s, k in zip(state, k3)], parameters)
    
    newState = [s + (stepSize / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i) 
                for s, k1_i, k2_i, k3_i, k4_i in zip(state, k1, k2, k3, k4)]
    
    return newState
