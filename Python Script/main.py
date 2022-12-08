# make a header comment for me to type in
""""
This is a script for part 2 of the fall 2022 dynamics project.
It will be used to calculate the minimum drop height for a ball to complete a loop on a ramp.

The true ramp angle is unknown, so we will just use 30 degrees for now.
The true radius of the loop is unknown, so we will just use 1 meter for now.
The mass of the ball is .028 kg
The radius of the ball .0159 m
The density of the air is 1.225 kg/m^3
The drag coefficient of a sphere is 0.47
The ball is rolling in between two rails, .011 m apart

We will use the bisection method to find the minimum drop height.
We will use a tolerance of 0.001 m.

We will use a timestep of 0.001 s.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

# define global constants
MASS = 0.028  # mass of the ball in kg
RADIUS = .5  # radius of the loop in m
RAMP_ANGLE = 50  # angle of the ramp in degrees
theta_rad = RAMP_ANGLE * math.pi / 180
theta_ball = 47.795
theta_ball_rad = theta_ball * math.pi / 180
g = 9.8  # gravitational acceleration in m/s^2
rho = 1.225  # density of the air in kg/m^3
Cd = 0.47  # drag coefficient of a sphere
r_ball = 0.0159  # radius of the ball in m
mu_s = .55


def bisection(f, a, b, tol):
    while abs(a - b) > tol:
        c = (a + b) / 2
        if f(c) == 0:
            return c
        elif f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2


def f(h):
    N = (MASS * g * math.cos(theta_rad)) / (2 * math.cos(theta_ball_rad))
    term1 = MASS * g * h
    term2 = (mu_s * N * (h - RADIUS * (1 - math.cos(theta_rad))) / (math.sin(theta_rad))) * (RADIUS * (theta_rad + math.pi))
    term3 = -MASS * g * 2 * RADIUS
    term4 = -(MASS * RADIUS * g) / 2
    term5 = (g * MASS * r_ball ** 2) / 5
    VALUE = term1 + term2 + term3 + term4 + term5
    return VALUE


def main():
    a = .001
    b = 100
    tol = .001

    min_drop_height = bisection(f, a, b, tol)

    print(f"The minimum drop height is {min_drop_height} m.")
    print(f"The minimum drop height is {min_drop_height * 39.3701} in.")


if __name__ == '__main__':
    main()
