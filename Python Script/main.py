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
RADIUS = 5 * 0.0254  # radius of the loop in inches converted to meters
RAMP_ANGLE = 50  # angle of the ramp in degrees
theta_rad = RAMP_ANGLE * math.pi / 180
g = 9.8  # gravitational acceleration in m/s^2


def bisection(f, a, b, tol, params):
    iters = 0
    while abs(a - b) > tol:
        iters += 1
        c = (a + b) / 2
        if f(c, params) == 0:
            return c
        elif f(a, params) * f(c, params) < 0:
            b = c
        else:
            a = c
    return (a + b) / 2, iters


def f(h, params):
    MASS = params[0]
    r_ball = params[1]
    theta_ball_rad = params[2]
    mu_s = params[3]

    N = (MASS * g * math.cos(theta_rad)) / (2 * math.cos(theta_ball_rad))
    term1 = MASS * g * h
    term2 = (mu_s * N * (h - RADIUS * (1 - math.cos(theta_rad))) / (math.sin(theta_rad)))
    term3 = mu_s * N * (RADIUS * (theta_rad + math.pi))
    term4 = MASS * g * 2 * RADIUS
    term5 = (MASS * RADIUS * g) / 2
    term6 = (g * MASS * r_ball ** 2) / (5 * RADIUS)
    VALUE = term1 - term2 - term3 - term4 - term5 - term6
    return VALUE


def main():
    a = .001
    b = 100
    tol = .0000001

    MASS = 0.028  # mass of the ball in kg
    r_ball = 0.0159  # radius of the ball in m
    theta_ball = 47.795
    theta_ball_rad = theta_ball * math.pi / 180
    mu_s = .55

    params = [MASS, r_ball, theta_ball_rad, mu_s]

    min_drop_height_rubber, iters = bisection(f, a, b, tol, params)

    print(f'\nRUBBER BALL')
    print(f'The radius of the loop is {RADIUS} m')
    print(f"The minimum drop height is {min_drop_height_rubber} m.")
    print(f"The minimum drop height is {min_drop_height_rubber * 39.3701} in.")
    print(f'The number of iterations is {iters}.')

    MASS = .003
    r_ball = .015
    theta_ball = 51.06
    theta_ball_rad = theta_ball * math.pi / 180
    mu_s = .5

    params = [MASS, r_ball, theta_ball_rad, mu_s]

    min_drop_height_plastic, iters = bisection(f, a, b, tol, params)

    print(f'\nPLASTIC BALL')
    print(f'The radius of the loop is {RADIUS} m')
    print(f"The minimum drop height is {min_drop_height_plastic} m.")
    print(f"The minimum drop height is {min_drop_height_plastic * 39.3701} in.")
    print(f"The number of iterations is {iters}.")

    MASS = .011
    r_ball = .014
    theta_ball = 55.44
    theta_ball_rad = theta_ball * math.pi / 180
    mu_s = .35

    params = [MASS, r_ball, theta_ball_rad, mu_s]

    min_drop_height_steel, iters = bisection(f, a, b, tol, params)

    print(f'\nSTEEL BALL')
    print(f'The radius of the loop is {RADIUS} m')
    print(f"The minimum drop height is {min_drop_height_steel} m.")
    print(f"The minimum drop height is {min_drop_height_steel * 39.3701} in.")
    print(f"The number of iterations is {iters}.")


if __name__ == '__main__':
    main()
