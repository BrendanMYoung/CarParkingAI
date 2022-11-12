import numpy as np
import scipy
import matplotlib.pyplot as plt
import math
from scipy import interpolate
population_size = 1
mutation_rate = .005
size_of_parameter_vector = 7
max_generations = 100
K = 200

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
def newChromeosome():
    return np.random.choice(chromosomes, size_of_parameter_vector).tolist()
# Creates a new sequence of chromosomes for an entire time step
def newChromesomoneSequence():
    chromeosomeSequence = []
    for i in range(0, total_entries, 1):
        chromeosomeSequence.append(newChromeosome())

    return chromeosomeSequence
# Generates the entire initial population percept sequences
def generateNewRandomPopulation():
    population = []
    for i in range(0, population_size, 1):
        population.append(newChromesomoneSequence())

    return population

# Converts the number from binary to decimal, the first index determines if the number is negative or positive
def binaryToDecimal(binary_array):
    decimal = 0
    for i in range(len(binary_array) -1, -1, -1):
        decimal = decimal + binary_array[i] * pow(2, len(binary_array)- (i + 1))
    if binary_array[0] == 0:
        decimal *= -1
    return decimal
# Converts a given binary input into standardized accelerations and headings
# Convert from - 2^(size_of_parameter_vector - 1), to [min, max]
max = pow(2, size_of_parameter_vector)
def convertHeading(binary_chromosome):
    oldVal = binaryToDecimal(binary_chromosome)
    newRange = max_heading - min_heading
    return (((oldVal - 0) * newRange ) / max) + min_heading
# Converts acceleration
def convertAcceleration(binary_chromosome):
    oldVal = binaryToDecimal(binary_chromosome)
    newRange = max_acceleration - min_acceleration
    return (((oldVal - 0) * newRange) / max) + min_acceleration

def convertChromesomeSequence(chromesome):
    acceleration_history = []
    heading_history = []
    x = []
    y = []
    for i in range(0, len(chromesome), 2):

        heading = convertHeading(chromesome[i])
        acceleration = convertAcceleration(chromesome[i+1])

        acceleration_history.append(heading)
        heading_history.append(acceleration)

    time = np.linspace(0, 10, num=10, endpoint = True)
    s = interpolate.splrep(time, acceleration_history, s=0)
    h = interpolate.splrep(time, heading_history, s=0)

    acceleration_new = interpolate.splev(all_time_steps_for_x, s, der=0)
    heading_new = interpolate.splev(all_time_steps_for_x, h, der=0)
    #tck = interpolate.splrep(time, acceleration_history, s=0)
    #print(tck)
    plt.figure()
    for timeStep in range(0, len(all_time_steps_for_x)):
        current_velocity = acceleration_new[timeStep]
        current_heading = heading_new[timeStep]

        x.append( current_heading* math.cos(current_velocity))
        y.append( current_heading * math.sin(current_velocity))
    plt.plot(x, y)
    plt.show()
# Basically converts an entire
def convertPopulationSeuqence(pops):
    for pop in pops:
        convertChromesomeSequence(pop)

def isInBounds(x, y):
    if x <= -4 and y > 3:
        return True
    elif x >-4 and x< 4 and y > -1:
        return True
    elif x>= 4 and y > 3:
        return True
    return False

def main():
    generateNewRandomPopulation()
    convertPopulationSeuqence(generateNewRandomPopulation())



if __name__ == "__main__":
    main()