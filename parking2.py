import numpy as np
import sys
import scipy
import random
import matplotlib.pyplot as plt
import math
from scipy import interpolate
population_size = 201
mutation_rate = .005
size_of_parameter_vector = 7
max_generations = 1200

min_heading = -0.524
max_heading = 0.524

min_acceleration = -5
max_acceleration = 5

time_steps = 10
steps_per_time_step = 0.1
# List of the time steps to use for interpolation
all_time_steps_for_x = np.arange(0, time_steps, steps_per_time_step)

amt_control_variables = 2
total_entries = time_steps * amt_control_variables
# Creates a new chromosome for a single time step
chromosomes = [0, 1]
total_chromosomes = 10

initial_state =  [0, 8, 0, 0]
final_state = np.array([0, 0, 0, 0])
K = 200
cost_for_being_out_of_bounds = 1
stop_criteria = 0.1

total_bits = size_of_parameter_vector

def newChromeosome():
    return np.random.choice(chromosomes, size_of_parameter_vector).tolist()
def createNewChromeosome(val):
    chromesome = [1, 0, 0, 0, 0, 0, 0]
    return chromesome
