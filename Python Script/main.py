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
RADIUS = 1.0  # radius of the loop in m
RAMP_ANGLE = 30.0  # angle of the ramp in degrees
g = 9.8  # gravitational acceleration in m/s^2
rho = 1.225  # density of the air in kg/m^3
Cd = 0.47  # drag coefficient of a sphere


def ball_motion(min_height, max_height, tolerance, dt):
    

def main():
    gaming = 0


if __name__ == '__main__':
    main()
