import numpy as np
import time
import random
import matplotlib.pyplot as plt
import math
from scipy import interpolate
population_size = 451
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

total_bits = size_of_parameter_vector * (total_entries - 2)  #subtract 2 for the initial 2 states

# Creates a new gray coded binary sequence
def newChromosome():
    chromosomeSequence = [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    chromosomeSequence += np.random.choice(chromosomes, total_bits).tolist()
    return chromosomeSequence

# Take in an array of 7 values for a value sequence of a gene to convert it to decimal
def binaryToDecimal(binary_array):
    decimal = 0
    for i in range(len(binary_array) -1, -1, -1):
        decimal = decimal + binary_array[i] * pow(2, len(binary_array)- (i + 1))
    if binary_array[0] == 0:
        decimal *= -1
    return decimal

def flip(c):
    return 1 if(c == 0) else 0;
def graytoBinary(gray):
    binary = [];

    binary.append(gray[0])

    # Compute remaining bits
    for i in range(1, len(gray)):

        # If current bit is 0,
        # concatenate previous bit
        if (gray[i] == 0):
            binary.append(binary[i - 1])

        # Else, concatenate invert
        # of previous bit
        else:
            binary.append(flip(binary[i - 1]))

    return binary;

# Converts a given binary input into standardized accelerations and headings
# Convert from - 2^(size_of_parameter_vector) - 1, to [min, max]
max = pow(2, size_of_parameter_vector)-1
def convertHeading(binary_chromosome):
    #binary_chromosome = graytoBinary(graycode_chromosome.tolist())
    oldVal = abs(binaryToDecimal(binary_chromosome))
    newRange = max_heading - min_heading
    return ((oldVal * newRange ) / max) + min_heading
# Converts acceleration
def convertAcceleration(binary_chromosome):
    #binary_chromosome = graytoBinary(graycode_chromosome.tolist())
    oldVal = abs(binaryToDecimal(binary_chromosome))
    newRange = max_acceleration - min_acceleration
    return ((oldVal * newRange) / max) + min_acceleration

def costFunction(local_final_state, cf):
    if cf == 0:
        cost = np.linalg.norm( local_final_state - final_state)
    else:
        cost = K + cf
    return cost
def generateNewPopulation():
    population = []
    for i in range(0, population_size, 1):
        population.append(newChromosome())
    return population

def convertChromesomeSequence(non_split_chromosome):
    chromesome = np.array_split(non_split_chromosome, total_entries)
    optimization_parameter = []
    acceleration_history = []
    heading_history = []
    x = []
    y = []
    for i in range(0, len(chromesome), 2):

        heading = convertHeading(chromesome[i])
        acceleration = convertAcceleration(chromesome[i+1])

        acceleration_history.append(acceleration)
        heading_history.append(heading)
        optimization_parameter.append(heading)
        optimization_parameter.append(acceleration)

    time = np.linspace(0, time_steps, num=time_steps, endpoint = True)
    s = interpolate.CubicSpline(time, acceleration_history, bc_type='natural')
    h = interpolate.CubicSpline(time, heading_history, bc_type='natural')
    x_new = np.linspace(0, 10, 100)
    acceleration_new = s(x_new)
    heading_new = h(x_new)

    current_x = initial_state[0]
    current_y = initial_state[1]
    #print(heading_new)
    local_cost = 0
    for timeStep in range(0, len(all_time_steps_for_x)):
        current_velocity = acceleration_new[timeStep]
        current_heading = heading_new[timeStep]
        current_x += steps_per_time_step * current_velocity * math.cos(current_heading)
        current_y += steps_per_time_step * current_velocity * math.sin(current_heading)

        if isNotInBounds(current_x, current_y):
            local_cost += math.pow(getYValue(current_x) - current_y, 2 ) * steps_per_time_step
        x.append( current_x)
        y.append( current_y)
    local_final_state = np.array([current_x, current_y, current_heading, current_velocity])
    cost = costFunction(local_final_state, local_cost)
    return {
        "cost": cost,
        "fitness": 1 / (cost+1),
        "x": x,
        "y": y,
        "acceleration": acceleration_new,
        "heading": heading_new,
        "optimization_vector": optimization_parameter
    }
    #return [cost, x, y]

def crossOver(pop1, pop2):
    # For every chromsome in each population
    cross_over_pivot = np.random.randint(0, len(pop1))
    second_pivot = np.random.randint(cross_over_pivot, len(pop1))
    new_child1 = pop2[0:cross_over_pivot] + pop1[cross_over_pivot:second_pivot] + pop2[second_pivot:]
    new_child2 = pop1[0:cross_over_pivot] + pop2[cross_over_pivot:second_pivot] + pop1[second_pivot:]
    new_c_1 = []
    new_c_2 = []
    # do mutation of children
    for i in range(len(new_child1)):
        new_c_1.insert(i, random.choices(
            [new_child1[i], (new_child1[i] + 1 % 2)],
            weights=[(1 - mutation_rate), mutation_rate]
        )[0]
                          )
        new_c_2.insert(i, random.choices(
            [new_child2[i], (new_child2[i] + 1 % 2)],
            weights=[1 - mutation_rate, mutation_rate]
        )[0])
    return [new_c_1, new_c_2]
# Do crossover breeding
def newPopulationCrossover(chance_array, pops):
    new_chromeosomes = []
    for i in range(math.floor(population_size/2)):

        new_indexs = np.random.choice(
            population_size,
            2,
            p=chance_array
        )
        pop1 = pops[new_indexs[0]]
        pop2 = pops[new_indexs[1]]
        new_sequences = crossOver(pop1, pop2)
        new_chromeosomes.append(new_sequences[0])
        new_chromeosomes.append(new_sequences[1])

    return new_chromeosomes

# Basically converts an entire
def doPopulationGeneration(pops, generation):
    fitness_of_pops = []
    total_combined_fitness = 0
    most_fit_pop_index = 0
    for i in range(len(pops)):
        pop = pops[i]
        values = convertChromesomeSequence(pop) # converts chromeosequence to cost function
        total_combined_fitness += values["fitness"] # get the cost index
        fitness_of_pops.insert(i, values)
        if values["fitness"] > fitness_of_pops[most_fit_pop_index]["fitness"]:
            most_fit_pop_index = i
    print("Generation " +str(generation) +" : J = " + str(fitness_of_pops[most_fit_pop_index]["cost"]) )


    if fitness_of_pops[most_fit_pop_index]["cost"] <= stop_criteria or generation >= max_generations:
        return None, fitness_of_pops[most_fit_pop_index]

    # Calculate fitness chance and setup random selection for crossover with pops
    chance_array = []
    for i in range(len(fitness_of_pops)):
        pop = fitness_of_pops[i]
        chance_array.insert(i, pop["fitness"] / total_combined_fitness)

    new_generation = newPopulationCrossover(chance_array, pops)
    #print(new_generation)
    new_generation.append(pops[most_fit_pop_index])
    return new_generation, fitness_of_pops[most_fit_pop_index]

def finalState(most_fit):
    print("")
    print("Final state values:")
    print("x_f = " + str(most_fit["x"][-1]))
    print("y_f = " + str(most_fit["y"][-1]))
    print("alpha_f = " + str(most_fit["heading"][-1]))
    print("v_f = " + str(most_fit["acceleration"][-1]))

    figure, axis = plt.subplots(5, 1, figsize=(4, 2 * 5), tight_layout= True)

    axis[0].set_xlim([-15, 15])
    axis[0].set_ylim([-10, 15])
    axis[0].plot([-15, -4], [3, 3], 'k-', lw=1)
    axis[0].plot([-4, -4], [-1, 3], 'k-', lw=1)
    axis[0].plot([-4, 4], [-1, -1], 'k-', lw=1)
    axis[0].plot([4, 4], [-1, 3], 'k-', lw=1)
    axis[0].plot([15, 4], [3, 3], 'k-', lw=1)
    axis[0].plot(most_fit["x"], most_fit["y"])
    axis[0].plot([-15, -4], [3, 3], 'k-', lw=1)
    axis[0].plot([-4, -4], [-1, 3], 'k-', lw=1)
    axis[0].plot([-4, 4], [-1, -1], 'k-', lw=1)
    axis[0].plot([4, 4], [-1, 3], 'k-', lw=1)
    axis[0].plot([15, 4], [3, 3], 'k-', lw=1)
    plt.setp(axis[0], xlabel="x (ft)", ylabel="y (ft)")

    axis[1].plot(all_time_steps_for_x, most_fit["acceleration"])
    plt.setp(axis[1], xlabel="Time (s)", ylabel="beta (ft/s^2)")

    axis[2].plot(all_time_steps_for_x, most_fit["heading"])
    plt.setp(axis[2], xlabel="Time (s)", ylabel="alpha (theta/s^2)")

    axis[3].plot(all_time_steps_for_x, most_fit["x"])
    plt.setp(axis[3], xlabel="Time (s)", ylabel="x (ft)")

    axis[4].plot(all_time_steps_for_x, most_fit["y"])
    plt.setp(axis[4], xlabel="Time (s)", ylabel="y (ft)")

    plt.show()

def isNotInBounds(x, y):
    if x <= -4 and y > 3:
        return False
    elif (x >-4 and x< 4) and y > -1:
        return False
    elif x>= 4 and y > 3:
        return False
    return True

def getYValue(x):
    if x <= -4:
        return 3
    elif x > -4 and x < 4:
        return -1
    elif x >= 4:
        return 3
    return 0

def main():
    start_time = time.time()
    generations = 0
    current_population = generateNewPopulation()
    while (current_population != None and time.time() -start_time < 420):
        current_population, most_fit = doPopulationGeneration(current_population, generations)
        generations += 1

    finalState(most_fit)
    file = open('control.dat', 'w')
    for val in most_fit["optimization_vector"]:
        file.write(str(val)+"\n")
    file.close()

if __name__ == "__main__":
    main()