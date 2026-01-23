from geodesic import geodesicEquations
from integrator import integrateGeodesic
import math
import matplotlib.pyplot as plt

SOLAR_MASS_KM = 1.476625
C_KM_S = 299_792.458

inputStr = input(
    "Enter values separated by spaces:\n"
    "coordinateSystem(g=geometric, r=real) "
    "M initialRadius initialAngle emissionAngleFraction "
    "maxAffineParameter initialStepSize outputInterval pathOutput(y/n)\n"
    "Example: g 1.0 6.0 0.0 0.5 50 0.01 50 y\n"
)

inputs = inputStr.strip().split()
if len(inputs) != 9:
    raise ValueError("Expected exactly 9 inputs.")

coordSystem = inputs[0].lower()
rawMass = float(inputs[1])
initialRadius = float(inputs[2])
initialAngle = float(inputs[3])
emissionAngleFraction = float(inputs[4])
maxAffineParameter = float(inputs[5])
initialStepSize = float(inputs[6])
outputInterval = int(inputs[7])
recordTrajectory = inputs[8].lower() == 'y'

if coordSystem == 'r':
    blackHoleMass = rawMass * SOLAR_MASS_KM
    initialAngle = math.radians(initialAngle)
    cFactor = C_KM_S
elif coordSystem == 'g':
    blackHoleMass = rawMass
    cFactor = 1.0
else:
    raise ValueError("coordinateSystem must be 'g' or 'r'.")

emissionAngle = emissionAngleFraction * math.pi
lapse = 1 - 2 * blackHoleMass / initialRadius

localTimeComponent = 1.0
localRadialComponent = math.cos(emissionAngle)
localAngularComponent = math.sin(emissionAngle)

timeVelocity = localTimeComponent / math.sqrt(lapse)
radialVelocity = localRadialComponent * math.sqrt(lapse)
angularVelocity = localAngularComponent / initialRadius

initialState = [
    0.0,
    initialRadius,
    initialAngle,
    timeVelocity,
    radialVelocity,
    angularVelocity
]

# --- Setup live plot ---
plt.ion()
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(-20, 20)
ax.set_ylim(-20, 20)
line, = ax.plot([], [], 'b-')
trailLength = 20  # number of points to keep in trail
xs, ys = [], []

def plotCallback(x, y):
    xs.append(x)
    ys.append(y)
    if len(xs) > trailLength:
        xs.pop(0)
        ys.pop(0)
    line.set_data(xs, ys)
    plt.draw()
    plt.pause(0.001)

trajectory, status, finalTime = integrateGeodesic(
    geodesicEquations,
    initialState,
    blackHoleMass,
    maxAffineParameter,
    initialStepSize,
    outputInterval,
    recordTrajectory,
    plotCallback=plotCallback
)

timeOfFlight = finalTime / cFactor * 1e9 if coordSystem == 'r' else finalTime

plt.ioff()
plt.show()

if recordTrajectory:
    print("Path output:")
    print(", ".join(f"({x:.5f}, {y:.5f})" for x, y in trajectory))

print(f"Status: {status}")
lastX, lastY = trajectory[-1] if trajectory else (
    initialRadius * math.cos(initialAngle),
    initialRadius * math.sin(initialAngle)
)
unitLabel = "km" if coordSystem == 'r' else "geometric units"
print(f"Last observed coordinates ({unitLabel}): ({lastX:.5f}, {lastY:.5f})")
timeLabel = "ns" if coordSystem == 'r' else "geometric units"
print(f"Time of flight ({timeLabel}): {timeOfFlight:.5f}")
