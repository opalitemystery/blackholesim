from geodesic import geodesicEquations
from integrator import integrateGeodesic
import math

SOLAR_MASS_KM = 1.476625  # 1 solar mass = 1.476625 km
C_KM_S = 299_792.458      # speed of light in km/s

inputStr = input(
    "Enter values separated by spaces:\n"
    "coordinateSystem(Modes: g=geometric, r=real) "
    "M initialRadius initialAngle initialRadialVelocity initialAngularVelocity "
    "maxAffineParameter initialStepSize outputInterval pathOutput(y/n)\n"
    "Example: r 1.0 3.01 0.0 0.0 0.192 50 0.01 50 y\n"
)

inputs = inputStr.strip().split()
if len(inputs) != 10:
    raise ValueError("Expected exactly 10 inputs.")

coordSystem = inputs[0].lower()
rawMass = float(inputs[1])
rawRadius = float(inputs[2])
rawAngle = float(inputs[3])
initialRadialVelocity = float(inputs[4])
initialAngularVelocity = float(inputs[5])
maxAffineParameter = float(inputs[6])
initialStepSize = float(inputs[7])
outputInterval = int(inputs[8])
recordTrajectory = inputs[9].lower() == 'y'

# --- Scale inputs according to chosen system ---
if coordSystem == 'r':  # real units
    blackHoleMass = rawMass * SOLAR_MASS_KM  # solar mass → km
    initialRadius = rawRadius
    initialAngle = math.radians(rawAngle)    # degrees → radians
    c_factor = C_KM_S                         # for time conversion
elif coordSystem == 'g':  # geometric units
    blackHoleMass = rawMass
    initialRadius = rawRadius
    initialAngle = rawAngle
    c_factor = 1.0
else:
    raise ValueError("coordinateSystem must be 'g' or 'r'.")

# --- Compute initial time velocity using null condition ---
lapse = 1 - 2*blackHoleMass/initialRadius
initialTimeVelocity = math.sqrt(
    (initialRadius**2 * initialAngularVelocity**2 + initialRadialVelocity**2 / lapse) / lapse
)

initialState = [
    0.0,                 # t
    initialRadius,       # r
    initialAngle,        # phi
    initialTimeVelocity, # tDot
    initialRadialVelocity, # rDot
    initialAngularVelocity # phiDot
]

# --- Run integration ---
trajectory, status, finalTime = integrateGeodesic(
    geodesicEquations,
    initialState,
    blackHoleMass,
    maxAffineParameter,
    initialStepSize,
    outputInterval,
    recordTrajectory
)

# --- Convert time to nanoseconds if real units ---
timeOfFlight = finalTime / c_factor * 1e9 if coordSystem == 'r' else finalTime

# --- Print results ---
if recordTrajectory:
    print("Path output:")
    print(", ".join(f"({x:.5f}, {y:.5f})" for x, y in trajectory))

print(f"Status: {status}")
lastX, lastY = trajectory[-1] if trajectory else (initialRadius*math.cos(initialAngle), initialRadius*math.sin(initialAngle))
unitLabel = "km" if coordSystem == 'r' else "geometric units"
print(f"Last observed coordinates ({unitLabel}): ({lastX:.5f}, {lastY:.5f})")
timeLabel = "ns" if coordSystem == 'r' else "geometric units"
print(f"Time of flight ({timeLabel}): {timeOfFlight:.5f}")
