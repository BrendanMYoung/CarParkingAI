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

initial_state =  [0, 8, 0, 0]
final_state = np.array([0, 0, 0, 0])
K = 200
cost_for_being_out_of_bounds = 1
stop_criteria = 0.1

def newChromeosome():
    return np.random.choice(chromosomes, size_of_parameter_vector).tolist()
def createNewChromeosome(val):
    chromesome = [1, 0, 0, 0, 0, 0, 0]
    return chromesome

# Creates a new sequence of chromosomes for an entire time step
def newChromesomoneSequence():
    chromeosomeSequence = []
    # Set up first chromeosome to have the initial state as timestep = 0
    chromeosomeSequence.insert(0, createNewChromeosome(initial_state[2]))
    chromeosomeSequence.insert(1, createNewChromeosome(initial_state[3]))
    # Generate random chromeosomes
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
# Convert from - 2^(size_of_parameter_vector) - 1, to [min, max]
max = pow(2, size_of_parameter_vector)-1
def convertHeading(binary_chromosome):
    oldVal = binaryToDecimal(binary_chromosome)
    newRange = max_heading - min_heading
    return (((oldVal - 0) * newRange ) / max) + min_heading
# Converts acceleration
def convertAcceleration(binary_chromosome):
    oldVal = binaryToDecimal(binary_chromosome)
    newRange = max_acceleration - min_acceleration
    return (((oldVal - 0) * newRange) / max) + min_acceleration

def costFunction(local_final_state, cf):
    if cf == 0:
        cost = np.linalg.norm( local_final_state - final_state)
    else:
        cost = K + cf
    return cost

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
    time = np.linspace(0, time_steps, num=time_steps+1, endpoint = True)
    #s = interpolate.splrep(time, acceleration_history, s=0)
   # h = interpolate.splrep(time, heading_history, s=0)
    #acceleration_new = interpolate.splev(all_time_steps_for_x, s, der=0)

    #heading_new = interpolate.splev(all_time_steps_for_x, h, der=0)
    s = interpolate.CubicSpline(time, acceleration_history, bc_type='natural')
    h = interpolate.CubicSpline(time, heading_history, bc_type='natural')
    x_new = np.linspace(0, 10, 100)
    acceleration_new = s(x_new)

    heading_new = h(x_new)

    #tck = interpolate.splrep(time, acceleration_history, s=0)
    #print(tck)

    current_x = initial_state[0]
    current_y = initial_state[1]
    #print(heading_new)
    local_cost = 0
    for timeStep in range(0, len(all_time_steps_for_x)):
        current_velocity = acceleration_new[timeStep]
        current_heading = heading_new[timeStep]
        current_x += steps_per_time_step * current_heading * math.cos(current_velocity)
        current_y += steps_per_time_step * current_heading * math.sin(current_velocity)
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
        "heading": heading_new
    }
    #return [cost, x, y]

def crossOver(pop1, pop2):
    # For every chromsome in each population
    cross_over_pivot = np.random.randint(0, len(pop1))
    new_child1 = pop2[0:cross_over_pivot] + pop1[cross_over_pivot:]
    new_child2 = pop1[0:cross_over_pivot] + pop2[cross_over_pivot:]
    new_c_1 = []
    new_c_2 = []
    # do mutation of children
    for i in range(len(new_child1)):
        new_c_1.insert(i, [])
        new_c_2.insert(i, [])

        for j in range(len(new_child1[i])):
            new_c_1[i].insert(j, random.choices(
                [new_child1[i][j], (new_child1[i][j]+1 % 2)],
                weights=[(1-mutation_rate), mutation_rate]
            )[0]
                              )
            new_c_2[i].insert(j, random.choices(
                [new_child2[i][j], (new_child2[i][j]+1 % 2)],
                weights=[1-mutation_rate, mutation_rate]
            )[0] )
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
    print("Generation " +str(generation) +" : " + str(fitness_of_pops[most_fit_pop_index]["cost"]) +" / " + str(fitness_of_pops[most_fit_pop_index]["fitness"]))


    if fitness_of_pops[most_fit_pop_index]["fitness"] >= 1-stop_criteria or generation >= max_generations:
        print('done')
        figure, axis = plt.subplots(3, 2)
        axis[0, 0].plot(fitness_of_pops[most_fit_pop_index]["x"], fitness_of_pops[most_fit_pop_index]["y"])
        axis[0, 0].set_title("x/y history")
        axis[0, 1].plot(all_time_steps_for_x, fitness_of_pops[most_fit_pop_index]["acceleration"])
        axis[0, 1].set_title("acceleration")
        axis[1, 0].plot(all_time_steps_for_x, fitness_of_pops[most_fit_pop_index]["heading"])
        axis[1, 0].set_title("heading")
        axis[2, 0].plot(all_time_steps_for_x, fitness_of_pops[most_fit_pop_index]["x"])
        axis[2, 0].set_title("x")
        axis[2, 1].plot(all_time_steps_for_x, fitness_of_pops[most_fit_pop_index]["y"])
        axis[2, 1].set_title("y")


        plt.show()
        return None

    # Calculate fitness chance and setup random selection for crossover with pops
    chance_array = []
    for i in range(len(fitness_of_pops)):
        pop = fitness_of_pops[i]
        chance_array.insert(i, pop["fitness"] / total_combined_fitness)

    new_generation = newPopulationCrossover(chance_array, pops)
    #print(new_generation)
    new_generation.append(pops[most_fit_pop_index])
    return new_generation

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
    generations = 0
    current_population = generateNewRandomPopulation()
    most_fit = None
    while (current_population != None):
        current_population = doPopulationGeneration(current_population, generations)
        generations +=1



if __name__ == "__main__":
    main()