import math

Gravity = Lightspeed = Mass = 1
Mass = 1

def schwarzchildFactor(radius):
    return 1 - (2*Mass)/radius

def angleRate(L ,radius):
    return L/radius**2

def timeRate(E, radius):
    return E / schwarzchildFactor(radius)