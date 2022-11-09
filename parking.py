import numpy as np
import scipy
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d

population_size = 201
mutation_rate = .005
size_of_parameter_vector = 7
max_generations = 100
K = 200

min_heading = -0.524
max_heading = 0.524

min_acceleration = -5
max_acceleration = 5

time_steps = 10

# Creates a new chromosome for a single time step
chromosomes = [0, 1]
def newChromeosome():
    return np.random.choice(chromosomes, size_of_parameter_vector).tolist()
# Creates a new sequence of chromosomes for an entire time step
def newChromesomoneSequence():
    chromeosomeSequence = []
    for i in range(0, time_steps, 1):
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
    print(chromesome)
    outputArray = []
    for i in range(0, len(chromesome), 2):
        outputArray.append(convertHeading(chromesome[i]))
        outputArray.append(convertAcceleration(chromesome[i+1]))
    print(outputArray)

# Basically converts an entire
def convertPopulationSeuqence(pops):
    for pop in pops:
        convertChromesomeSequence(pop)

def main():
    generateNewRandomPopulation()
    convertPopulationSeuqence(generateNewRandomPopulation())
if __name__ == "__main__":
    main()