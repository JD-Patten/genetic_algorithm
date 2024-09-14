
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
import warnings
from scipy.spatial.transform import Rotation as R


class Parameter():

    def __init__(self,  name, limits = None, value = None):

        self.name = name
        self.limits = limits 
        self.value = value
        self.genotype_range = None

    def __repr__(self):
        return f"Parameter('{self.name}', {self.limits}, {self.value})"
    
    def validate(self):

        #limits and value can't both be None. 
        #If they are there's no way to define its genes or use it in any way
        if self.limits is None and self.value is None:
            raise ValueError(f"no limits or value given for parameter: {self.name}")
        
        # if limits are not provided, it's assumed that the parameter will not be used in the genotype
        # it will not be represented in binary and not used in any of the GA operations.
        if self.limits == None:
            self.used_for_genotype = False
        else:
            self.used_for_genotype = True

class Substructure():
    def __init__(self, type, name, limits_dict):
        self.type = type.lower()
        self.name = name
        self.limits_dict = limits_dict
        self.parameters: List[Parameter] = []

        if type.lower() == "dof":

            # currently not doing anything special for dof type substructures
            for param_name, limits in limits_dict.items():
                self.parameters.append(Parameter(limits = limits, name = param_name))

        else:
            for param_name, limits in limits_dict.items():
                self.parameters.append(Parameter(limits = limits, name = param_name))
                
        
    def __repr__(self):
        return f"Substructure('{self.type}', '{self.name}', {self.limits_dict})"
    
    def solve_sin_function(self, time):

        if self.type != 'dof':
            return None

        # Create a dictionary for parameter values
        param_dict = {param.name: param.value for param in self.parameters}

        # Extract parameters with default values if not provided
        amplitude = param_dict.get('amplitude')
        period = param_dict.get('period')
        time_shift = param_dict.get('time_shift')
        position_shift = param_dict.get('position_shift')


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

        # check to make sure each parameter is valid and give each parameter a range  where its genes 
        # will be found in the genotype

        count = 0
        for substructure in self.substructures:
            for parameter in substructure.parameters:
                parameter.validate()                                        
                if parameter.used_for_genotype:
                    parameter.genotype_range = (count, count + resolution)
                    count += resolution
                else: 
                    parameter.genotype_range = None

        # if no specific genotype is given, create randomized genotype from the substructures
        if genotype == None:
            self.randomize_genes()
        else:
            self.genotype = genotype

        self.update_phenotype()

        self.set_transforms()
    
    def __repr__(self):
        return f"Organism({self.substructures}, {self.resolution}, {self.genotype})"
    
    def update_phenotype(self):
        """
        updates the organisms' parameter values based on the current genotype
        """
    
        resolution = self.resolution

        for substructure in self.substructures:
            for parameter in substructure.parameters:
                if parameter.used_for_genotype:
    
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
                if parameter.used_for_genotype:

                    # Set random genes for every parameter that is being used in the genotype
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

    def solve_inverse_kinematics(self, time):

        l1 = self.l1
        l2 = self.l2

        #get the end effector pose based on organism's sin waves
        end_effector_pose = self.solve_sin_functions(time)

        translation = [end_effector_pose['x'] * 0.001,
                       end_effector_pose['y'] * 0.001, 
                       end_effector_pose['z'] * 0.001]
        
        rotation = [end_effector_pose['roll'], 
                    end_effector_pose['pitch'], 
                    end_effector_pose['yaw']]
        
        rotation_matrix = R.from_euler('xyz', rotation, degrees = True)

        #solve the inverse kinematics for the 6 arms
        angles = []
        for arm_number in range(6):

            #rotate the end connection point
            end_connection_point = rotation_matrix.apply(self.end_effector_arm_connections[arm_number])

            #translate the end connection point so it's at the correct position in the base frame
            end_connection_point += translation

            #rotate and translate the point so it is in the servo's frame
            end_connection_point -= self.servo_translations[arm_number]
            end_connection_point = self.servo_rotations[arm_number].apply(end_connection_point)

            x = end_connection_point[0]
            y = end_connection_point[1]
            z = end_connection_point[2]

            #adding and subtracting the acos portion decides which of the two possible solutions to choose
            #using if z < 0 uses different solutions for the upside down servos. this keeps the arms pointing outwards
            if z < 0:
                angle = math.atan2(z,y) + math.acos((l2**2 - y**2 - z**2 - x**2 - l1**2) / (-2 * l1 * math.sqrt(y**2 + z**2)))
            else:
                angle = math.atan2(z,y) - math.acos((l2**2 - y**2 - z**2 - x**2 - l1**2) / (-2 * l1 * math.sqrt(y**2 + z**2)))
        
            angles.append(angle)

        return angles

    def vary_speed(self, multiplier):
        # iterate through all the parameters varying them if they are a period parameter
        for substructure in self.substructures:
            for parameter in substructure.parameters:
                if parameter.name == "period":
                   
                    target_value = parameter.value * multiplier

                    if parameter.used_for_genotype: 
                        
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

                    # if the parameter is not being used in the genotype (probably because constant periods are being used)
                    # we just update the prameter value to the target value
                    else:
                        parameter.value *= multiplier


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
                    
                    parameter.value = period

                    if parameter.limits != None:
                        warnings.warn(f'setting a constant period on parameter ({parameter.name}) with limits ({parameter.limits})')

    def set_transforms(self):

        # Initialize an empty dictionary to store parameters
        ik_parameters_dict = {}

        # Iterate through Substructure to collect parameters
        for ss in self.substructures:
            if ss.name == 'ik':
                for parameter in ss.parameters:
                    # Add parameter name and value to the dictionary
                    ik_parameters_dict[parameter.name] = parameter.value

        # set the arm lengths
        self.l1 = ik_parameters_dict['l1']
        self.l2 = ik_parameters_dict['l2']

        # set the translations for the servos and end effector connection points for inverse kinematics
        self.servo_translations = find_servo_or_connection_point_coords(ik_parameters_dict['servo_offset_1'],
                                                                          ik_parameters_dict['servo_offset_2'],
                                                                          ik_parameters_dict['servo_height'])
        
        
        self.end_effector_arm_connections = find_servo_or_connection_point_coords(ik_parameters_dict['end_connection_offset_1'],
                                                                                  ik_parameters_dict['end_connection_offset_2'],
                                                                                  ik_parameters_dict['end_connection_z_offset'])


        #the rotations are set up 'backwards' with yzx and negative signs because we are going from the servo to the base frame
        draft_angle = ik_parameters_dict['draft_angle']
        self.servo_rotations = [R.from_euler('zyx', [-120,  draft_angle,    0], degrees=True),
                                R.from_euler('zyx', [-120,  draft_angle, -180], degrees=True),
                                R.from_euler('zyx', [   0,  draft_angle,    0], degrees=True),
                                R.from_euler('zyx', [   0,  draft_angle, -180], degrees=True),
                                R.from_euler('zyx', [ 120,  draft_angle,    0], degrees=True),
                                R.from_euler('zyx', [ 120,  draft_angle, -180], degrees=True)]

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

def find_servo_or_connection_point_coords(offset_1, offset_2, z_height):

    d1 = 1/6 * (offset_1 - offset_2)
    d2 = 1/3 * (offset_1 - offset_2)
    c = offset_1 - (2*d1)
    radius = math.sqrt(c**2 + (d2)**2 - 2*c*d2 * math.cos(math.radians(60)))

    rotate_120 = R.from_euler('z', 120, degrees=True)
    rotate_240 = R.from_euler('z', 240, degrees=True)

    number1_translations = [-1 * math.sqrt((radius ** 2) - ((offset_1/2) **2)), offset_1/2, z_height]
    number6_translations = [number1_translations[0], -1 * number1_translations[1], z_height]
    number2_translations = rotate_240.apply(number6_translations)
    number3_translations = rotate_240.apply(number1_translations)
    number4_translations = rotate_120.apply(number6_translations)
    number5_translations = rotate_120.apply(number1_translations)


    translations = [number1_translations, number2_translations, number3_translations,
                    number4_translations, number5_translations, number6_translations]
    
    return translations   

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





