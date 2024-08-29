# definitions from Goldberg GA book

# Organism - a possible solution being tested. aka 
# Genotype - total genetic package of an organism aka Structure
# Genes - individual characters (usually 1 or 0) that make up the Genotype
# Coding - How the genotype mapped to the phenotype
# Parameter - the real value that several genes convert into to create
# Phenotype - the list of parameters that a genotype converts to

# substructure - a group of parameters that tend to go together. 
# resolution - number of bits (or genes) per parameter. ex: a resolution of 6 means 2^6 possibilites for that parameter

import random
import math
from typing import List
import copy



class Parameter():


    def __init__(self,  name, limits, value = None):
        self.name = name
        self.limits = limits 
        self.value = value
        self.genotype_range = None
    
    def __repr__(self):
        return f"Parameter('{self.name}', {self.limits}, {self.value})"

class Substructure():
    def __init__(self, type, name, list_of_limits, custom_parameter_names = None):
        self.type = type.lower()
        self.name = name
        self.list_of_limits = list_of_limits
        self.parameters: List[Parameter] = []

        if type.lower() == "dof":

            self.parameters.append(Parameter(limits = list_of_limits[0], name = "amplitude"))
            self.parameters.append(Parameter(limits = list_of_limits[1], name = "period"))
            self.parameters.append(Parameter(limits = list_of_limits[2], name = "position shift"))
            self.parameters.append(Parameter(limits = list_of_limits[3], name = "time shift"))

        else:
            for i, name in enumerate(custom_parameter_names):
                self.parameters.append(Parameter(limits = list_of_limits[i], name = name))
                
            return
        
    def __repr__(self):
        return f"Substructure('{self.type}', '{self.name}', {self.list_of_limits})"
    
    def solve_sin_function(self, time):

        if self.type != 'dof':
            return None

        # Create a dictionary for parameter values
        param_dict = {param.name: param.value for param in self.parameters}

        # Extract parameters with default values if not provided
        amplitude = param_dict.get('amplitude')
        period = param_dict.get('period')
        time_shift = param_dict.get('time shift')
        position_shift = param_dict.get('position shift')


        # Calculate the value
        value = amplitude * math.sin((math.pi * 2 / period) * (time + time_shift * period)) + position_shift

        return value

class Organism():
    def __init__(self, substructures: List[Substructure] , resolution = 5, genotype = None, number = None):
 
        self.fitness = 0.0
        self.scaled_fitness = None

        self.substructures = substructures
        self.resolution = resolution

        self.number = number

        # find the range in the genotype where each parameter can be found
        # Using muyltiparameter, mapped, fixed-point coding as described in Goldber p.82
        count = 0
        for substructure in self.substructures:
            for parameter in substructure.parameters:
                parameter.genotype_range = (count, count + resolution)
                count += resolution

        # if no specific genotype is given, create randomized genotype from the substructures
        if genotype == None:
            self.randomize_genes()
        else:
            self.genotype = genotype

        self.update_phenotype()
    
    def __repr__(self):
        return f"Organism({self.substructures}, {self.resolution}, {self.genotype})"
    
    def update_phenotype(self):
        """
        updates the parameter values based on the current genotype
        """

        resolution = self.resolution

        for substructure in self.substructures:
            for parameter in substructure.parameters:

                #find the right binary string
                range = parameter.genotype_range
                genotype_section = self.genotype[range[0]:range[1]]

                #convert the binary to base 10 and map it to the parameters' range
                binary_value = "".join(map(str, genotype_section))
                base_10_value = int(binary_value, base=2)
                mapped_value =  base_10_value / (2 ** resolution - 1) * (parameter.limits[1] - parameter.limits[0]) + parameter.limits[0]
                parameter.value = mapped_value

    def mutate(self, mutation_frequency = None):
        # set default mutation frequency
        if mutation_frequency == None:
            mutation_frequency = 0.01
        
        # flip the gene if mutation probability is met
        for i, gene in enumerate(self.genotype):
            rand = random.random()
            if rand < mutation_frequency:
                if gene == 0:
                    self.genotype[i] = 1
                elif gene == 1:
                    self.genotype[i] = 0
                else:
                    print(f"gene not mutated due to non binary value:{gene}")

        # update the parameter values so they match the genotype
        self.update_phenotype()

        return
    
    def scale_fitness(self, a, b):
        self.scaled_fitness = a * self.fitness + b
    
    def randomize_genes(self):
        resolution = self.resolution

        self.genotype = []
        for substructure in self.substructures:
            for parameter in substructure.parameters:
                start_of_range = len(self.genotype)
                end_of_range = start_of_range + resolution

                for i in range(resolution):
                    self.genotype.append(random.choice((0,1)))

                parameter.genotype_range = (start_of_range, end_of_range)

        self.update_phenotype()
    
    def solve_sin_functions(self, time):
        
        solutions_dict = {}
        for substructure in self.substructures:
            if substructure.type == 'dof':
                solutions_dict[substructure.name] = substructure.solve_sin_function(time)

        return solutions_dict

    def vary_speed(self, multiplier):

        for substructure in self.substructures:
            for parameter in substructure.parameters:
                if parameter.name == "period":
                    
                    target_value = parameter.value * multiplier

                    # make sure target value is in range
                    if target_value < parameter.limits[0]:
                        target_value = parameter.limits[0]
                    elif target_value > parameter.limits[1]:
                        target_value = parameter.limits[1]

                    # map the adjusted value to the binary range (rounding up to the nearest int)
                    max_binary_mapped = 2 ** self.resolution - 1
                    percent = (target_value - parameter.limits[0]) / (parameter.limits[1] - parameter.limits[0])
                    mapped_value = math.ceil(percent * max_binary_mapped)
                    
                    #convert to binary
                    binary_value = bin(mapped_value)[2:]

                    #add zeros to the front of the binary value if it's too short
                    for i in range(self.resolution - len(binary_value)):
                        binary_value = '0'+binary_value

                    for i in range(parameter.genotype_range[0], parameter.genotype_range[1]):
                        gene = int(binary_value[i - parameter.genotype_range[0]])
                        self.genotype[i] = gene

        self.update_phenotype()
    
    def get_phenotype_dict(self):
        ss_dict = {}
        for substructure in self.substructures:
            param_dict = {}
            for parameter in substructure.parameters:
                param_dict[parameter.name] = parameter.value
            ss_dict[substructure.name] = param_dict

        return ss_dict

    def set_constant_period(self, period):

        for substructure in self.substructures:
            for parameter in substructure.parameters:
                if parameter.name == "period":
                    
                    target_value = period

                    # make sure target value is in range
                    if target_value < parameter.limits[0]:
                        target_value = parameter.limits[0]
                    elif target_value > parameter.limits[1]:
                        target_value = parameter.limits[1]

                    # map the adjusted value to the binary range (rounding up to the nearest int)
                    max_binary_mapped = 2 ** self.resolution - 1
                    percent = (target_value - parameter.limits[0]) / (parameter.limits[1] - parameter.limits[0])
                    mapped_value = math.ceil(percent * max_binary_mapped)
                    
                    #convert to binary
                    binary_value = bin(mapped_value)[2:]

                    #add zeros to the front of the binary value if it's too short
                    for i in range(self.resolution - len(binary_value)):
                        binary_value = '0'+binary_value

                    for i in range(parameter.genotype_range[0], parameter.genotype_range[1]):
                        gene = int(binary_value[i - parameter.genotype_range[0]])
                        self.genotype[i] = gene

        self.update_phenotype()



class Population():

    
    def __init__(self, organisms: List[Organism] = None, number = None):
        
        if organisms == None:
            self.organisms = []
        else:
            self.organisms = organisms
        
        self.number = number

    def __repr__(self):
        return f"Population({self.organisms}, {self.number})"

    def diversity(self):
        """
        Calculates the population diversity by finding the total Hamming Distance between all unique pairs of organisms

        Returns the total distance / the expected total distance for randomized genes

        Completely randomized populations should be close to 1

        Completely uniform populations will return 0
        """

        num_organisms = len(self.organisms)

        # Number of unique pairs
        num_pairs = num_organisms * (num_organisms - 1) / 2

        # Expected Hamming distance for each pair
        expected_hamming_distance = len(self.organisms[0].genotype) * 0.5

        # Total expected diversity score for randomized population
        expected_random_distance = num_pairs * expected_hamming_distance
        

        distance = 0
        
        # iterates through all unique pairs of organisms
        for i in range(num_organisms):
            for j in range(i + 1, num_organisms):
                org1 = self.organisms[i]
                org2 = self.organisms[j]

                # find how many genes differ between each org
                for a, b in zip(org1.genotype, org2.genotype):
                    if a != b:
                        distance += 1

        diversity = distance / expected_random_distance

        return diversity




    def total_fitness(self):
        """
        returns the sum of the individual organisms fitnesses
        """
        total_fitness = 0.0

        for org in self.organisms:
            total_fitness += org.fitness
        
        return total_fitness


    def mutate(population, mutation_frequency = None):
        
        for organism in population.organisms:
            organism.mutate(mutation_frequency)
    
        return 


def crossover(org1: Organism, org2: Organism, cut_point = None):

    # if no crossover point is given a random one will be selected
    if cut_point == None:
        cut_point = random.randint(0,len(org1.genotype))

    elif cut_point < 0 or cut_point > len(org1.genotype):
        raise Exception(f"cut point '{cut_point}' is out of allowable range: 0-{len(org1.genotype)} ")


    genotype1 = org1.genotype[:cut_point] + org2.genotype[cut_point:]
    genotype2 = org2.genotype[:cut_point] + org1.genotype[cut_point:]

    off1 =  Organism(substructures= org2.substructures, resolution=org2.resolution, genotype= genotype1)
    off2 =  Organism(substructures= org1.substructures, resolution=org1.resolution, genotype= genotype2)


    return off1, off2

def selection(population: Population, use_fitness_scaling: bool)->List[Organism]:

    # (Goldberg p.76) 
    # scale the fitness with linear scaling f' = af + b  
    # scale to keep the average fitness the same and the max a proportion of the average (f'_max = C_mult * f_average)
    c_mult = 1.75

    #find the fitness stats (min, max, average, total sum) for the population
    total_fitness = 0
    max_fitness = population.organisms[0].fitness
    min_fitness = population.organisms[0].fitness

    for organism in population.organisms:

        if organism.fitness < min_fitness:
            min_fitness = organism.fitness
        
        if organism.fitness > max_fitness:  
            max_fitness = organism.fitness
        
        total_fitness += organism.fitness   

    average_fitness = total_fitness / len(population.organisms) 


    # find the linear scaling coeffiecients a and b
    if min_fitness > (c_mult * average_fitness - max_fitness) / (c_mult - 1):   # test to see if normal scaling will send min fitness value negative
        # normal scaling
        delta = max_fitness - average_fitness
        a = (average_fitness * c_mult - average_fitness) / (delta)              # a is the 'slope' aka df'/ df
        b = average_fitness * (max_fitness - c_mult * average_fitness) / delta  # b is the 'y intercept' aka f' where f = 0
    else:
        # scale as much as possible
        delta = average_fitness - min_fitness
        a = average_fitness / delta
        b = -1 * min_fitness * average_fitness / delta

    for organism in population.organisms:
        organism.scale_fitness(a,b)
    
    # if we're not using fitness scaling make the scaled fitness equal the normal fitness
    # to do - make this function not do all the scaling just to change it if we're not using it
    if use_fitness_scaling == False:
        for organism in population.organisms:
            organism.scaled_fitness = organism.fitness


    # find the total fitness of the population so we can create
    # probabilities of selection based on relative fitness
    total_scaled_fitness = 0
    for organism in population.organisms:
        total_scaled_fitness += organism.scaled_fitness



    #start the selection 
    selections = []
    for i in range(len(population.organisms)):

        # select a random number in the total scaled fitness
        random_number = random.random() * total_scaled_fitness

        #spin the "roulette wheel" by adding up the fitnesses until the random number has been reached
        sum = 0
        for organism in population.organisms:
            sum += organism.scaled_fitness
            if sum > random_number:
                selections.append(organism)
                break
    
    return selections


if __name__ == '__main__':

                            #amplitude,   period,   position shift,   time shift
    x_limits =                [[0,30],     [2.0,10.0],     [-25,25],       [0,1]] 
    y_limits =                [[0,30],     [2.0,10.0],     [-25,25],       [0,1]]  
    z_limits =                [[0,30],     [2.0,10.0],     [180,230],      [0,1]] 
    roll_limits =             [[0,20],     [2.0,10.0],     [-15,15],       [0,1]]
    pitch_limits =            [[0,20],     [2.0,10.0],     [-15,15],       [0,1]]  
    yaw_limits =              [[0,20],     [2.0,10.0],     [-15,15],       [0,1]]


    
    orgs = []
    for i in range(25):

        substructures = [Substructure('dof', 'x', x_limits),
                        Substructure('dof', 'y', y_limits),
                        Substructure('dof', 'z', z_limits),
                        Substructure('dof', 'roll', roll_limits),
                        Substructure('dof', 'pitch', pitch_limits),
                        Substructure('dof', 'yaw', yaw_limits)]
        
        orgs.append(Organism(substructures, resolution= 4))
 
    p = Population(orgs, 4)




    #print(f"genotype: {p.organisms[1].genotype}")
    for org in p.organisms:
        print(org.genotype)





