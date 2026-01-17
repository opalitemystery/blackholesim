import math

from integrator import rk4Step
from geodesic import geodesicRHS
from schwarzchild import Mass




deltaLambda = 0.01        
maxSteps = 20000          
captureRadius = 2 * Mass  
recordStep = 5            




state = [
    0.0,   
    10.0,  
    0.0    
]


parameters = {
    "E": 1.0,  
    "L": 5.0   
}



trajectory = []
lambdaValue = 0.0
status = "ESCAPED"

for step in range(maxSteps):

    t, r, angle = state


    if step % recordStep == 0:
        x = r * math.cos(angle)
        y = r * math.sin(angle)
        trajectory.append((x, y))

    if r <= captureRadius:
        status = "CAPTURED"
        break


    state = rk4Step(
        geodesicRHS,
        lambdaValue,
        state,
        deltaLambda,
        parameters
    )

    lambdaValue += deltaLambda




print("\nSimulation finished")
print("-------------------")
print(f"Status: {status}")
print(f"Final affine parameter λ: {lambdaValue:.4f}")

t, r, angle = state
x = r * math.cos(angle)
y = r * math.sin(angle)

print(f"Final position (r, angle): ({r:.6f}, {angle:.6f})")
print(f"Final position (x, y): ({x:.6f}, {y:.6f})")


trajectoryStr = "[" + ",".join(f"({x},{y})" for x, y in trajectory) + "]"

print("\nPhoton trajectory:")
print(trajectoryStr)

