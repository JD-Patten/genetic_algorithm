#!/usr/bin/env python3
import math
import random
import numpy as np
import rclpy
import itertools
import copy
 
from genetic_algorithm.msg import PopulationStats as PopulationStatsMsg, Organism as OrganismMsg, OrganismParameter as OrganismParameterMsg

from rclpy.node import Node   
from scipy.spatial.transform import Rotation as R

from geometry_msgs.msg import Pose, Point, Quaternion, PoseStamped
from sensor_msgs.msg import JointState
from std_msgs.msg import Float64, String
from rosgraph_msgs.msg import Clock

from ga_classes import *



servo_offset_1 = 95 * 0.001                   # mm to m Distance between servos 2-3, 4-5, 6-1
servo_offset_2 = 72 * 0.001                   # mm to m Distance between servos 1-2, 3-4, 5-6
servo_height = 43 * 0.001                     # mm to m 
    
end_connection_offset_1 = 50 * 0.001          # mm to m Distance between arm connection points 2-3, 4-5, 6-1
end_connection_offset_2 = 30 * 0.001          # mm to m Distance between arm connection points 1-2, 3-4, 5-6
end_connection_z_offset = -12.5 * 0.001       # mm to m 
    
draft_angle = 20                              # degrees
l1 = 65 * 0.001                               # mm to m 
    
                         #amplitude,   period,     position shift,   time shift
x_limits =                [[0,30],     [1.0,5.0],     [-25,25],       [0,1]] 
y_limits =                [[0,30],     [1.0,5.0],     [-25,25],       [0,1]]  
z_limits =                [[0,30],     [1.0,5.0],     [180,210],      [0,1]] 
roll_limits =             [[0,20],     [1.0,5.0],     [-15,15],       [0,1]]
pitch_limits =            [[0,20],     [1.0,5.0],     [-15,15],       [0,1]]  
yaw_limits =              [[0,20],     [1.0,5.0],     [-15,15],       [0,1]]


ik_genes_limits = [[40,100], [20,50], [150,250]]

planted_population = None

class SimulationManager(Node):
    def __init__(self):
        super().__init__("simulation_manager")

        # GA features
        self.use_planted_population = False
        self.use_fitness_scaling = False
        self.use_ik_genes = False                         # this toggles using additional genes to control the dimensions of the end effector and long arms
        self.use_varied_speed = False
        self.use_constant_period = False
        self.constant_period = 1.0
        self.check_for_hit_ground = False
        

        # GA tuning
        self.population_size = 6                         #number of organisms in each generation
        self.mutation_frequency = 0.01                    #probability a gene is flipped during mutation
        self.crossover_frequency = 0.85
        self.max_lifetime = 12                            #seconds
        self.ramp_up_time = 0.80                          #time for new organism to go from default parameters to it's own parameters
        self.resolution = 6                               #number of genes (1s or 0s) per parameter


        self.l2 =  194 *0.001                             # arm length in mm

        self.publish_joint_goal_states = False            #used to publish the arm angle goals being sent to the sim

        self.running = False
        self.sub = self.create_subscription(Clock,"/clock", self.get_sim_time, 1) 
        self.sub = self.create_subscription(PoseStamped, "/model/skipper/pose", self.get_robot_pose_from_sim, 1)
        self.sub = self.create_subscription(JointState, "/world/bot_world/model/skipper/joint_state", self.get_joint_states,1)
        self.angles_from_sim = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        self.distance_traveled_pub = self.create_publisher(Float64, "/distance_traveled", 1)
        self.population_info_pub = self.create_publisher(PopulationStatsMsg, "/pop_info",1)
        self.lifetime_pub = self.create_publisher(Clock,"/lifetime",1)
        self.bot_pose_pub = self.create_publisher(String, "/bot_pose",1)

        # publishers for the arm joints in the sim
        self.pub1 = self.create_publisher(Float64, "/model/skipper/joint/arm1_joint/_0/cmd_pos", 1)  
        self.pub2 = self.create_publisher(Float64, "/model/skipper/joint/arm2_joint/_0/cmd_pos", 1)  
        self.pub3 = self.create_publisher(Float64, "/model/skipper/joint/arm3_joint/_0/cmd_pos", 1)  
        self.pub4 = self.create_publisher(Float64, "/model/skipper/joint/arm4_joint/_0/cmd_pos", 1)  
        self.pub5 = self.create_publisher(Float64, "/model/skipper/joint/arm5_joint/_0/cmd_pos", 1)  
        self.pub6 = self.create_publisher(Float64, "/model/skipper/joint/arm6_joint/_0/cmd_pos", 1) 

        self.joint_goal_pub = self.create_publisher(JointState,"/arm_control_input", 1)

        self.timer = self.create_timer(0.03, self.set_end_effector_pose)

        self.sim_time = 0.0
        self.sim_time_offset = 0.0

        self.initialize_transforms()
        self.organism_starting_position = [0.0, 0.0, 0.0]
        self.create_initial_population()

        self.running = True

    def get_joint_states(self, msg):
        self.angles_from_sim = msg.position

    def publish_pop_info(self, population = None):

        if population == None:
            population = self.population


        msg = PopulationStatsMsg()
        msg.current_organism = self.current_organism_number
        msg.generation_number = self.current_generation_number
        msg.population_size = len(population.organisms)
        msg.population_fitness = population.total_fitness()

        msg.population_diversity = population.diversity()

        for organism in population.organisms:
            organism_msg = OrganismMsg()  
            organism_msg.max_distance = organism.fitness

            phenotype_dict = organism.get_phenotype_dict()
            print(phenotype_dict)
                
            dof_msg_dict = {'x': organism_msg.x, 'y': organism_msg.y, 'z': organism_msg.z,
                            'roll': organism_msg.roll, 'pitch': organism_msg.pitch, 'yaw': organism_msg.yaw}
            
            for dof, param_dict in phenotype_dict.items():
                
                param_msg = OrganismParameterMsg()

                param_msg.amplitude       = float(param_dict.get("amplitude"))
                param_msg.period          = float(param_dict.get("period"))
                param_msg.position_offset = float(param_dict.get("position shift"))
                param_msg.time_offset     = float(param_dict.get("time shift"))


                org_msg_dof = dof_msg_dict.get(dof)
                org_msg_dof.append(param_msg)


            msg.organisms.append(organism_msg)

        print(msg)
        self.population_info_pub.publish(msg)

    def initialize_transforms(self):
    
        # set the translations for the servos and end effector connection points for inverse kinematics
        self.servo_translations = self.find_servo_or_connection_point_coords(servo_offset_1, servo_offset_2, servo_height)
        self.end_effector_arm_connections = self.find_servo_or_connection_point_coords(end_connection_offset_1, end_connection_offset_2, end_connection_z_offset)

        #the rotations are set up 'backwards' with yzx and negative signs because we are going from the servo to the base frame
        self.servo_rotations = [R.from_euler('zyx', [-120,  draft_angle,    0], degrees=True),
                                R.from_euler('zyx', [-120,  draft_angle, -180], degrees=True),
                                R.from_euler('zyx', [   0,  draft_angle,    0], degrees=True),
                                R.from_euler('zyx', [   0,  draft_angle, -180], degrees=True),
                                R.from_euler('zyx', [ 120,  draft_angle,    0], degrees=True),
                                R.from_euler('zyx', [ 120,  draft_angle, -180], degrees=True)]
    
    def create_initial_population(self):
        
        if self.use_planted_population:
            self.population = planted_population
            self.population.number = 0
        
        else:

            organisms = []
            for i in range(self.population_size):

                # a new list of substructures is made for each new organism
                substructures = [Substructure('dof', 'x', x_limits),
                                Substructure('dof', 'y', y_limits),
                                Substructure('dof', 'z', z_limits),
                                Substructure('dof', 'roll', roll_limits),
                                Substructure('dof', 'pitch', pitch_limits),
                                Substructure('dof', 'yaw', yaw_limits)]

                validated = False
                while not validated:
                    org = Organism(substructures, resolution= self.resolution, number=i)
                    validated = self.check_parameters(org)

                organisms.append(org)

            self.population = Population(organisms, 0)

        self.get_logger().info(f'Initial Population: {self.population}') 
        
        #initialize all the values to track this population and future populations
        self.current_generation_number = 0
        self.current_organism_number = 0
        self.at_starting_position = True
        self.organism_position = [0, 0, 0]

        self.publish_pop_info(self.population)

    def check_parameters(self, organism: Organism):

        """
        Used to check if a max allowable angle and a max allowable
        angular velocity are respected by the parameters.

        Returns True if they are and False if they are not 

        If self.use_vary_speed is True, periods will be automatically 
        updated to match the max allowable angular velocity

        If self.use_constant_period is True, all periods will be reset
        before the check starts
        """


        hit_ground = False
        ground_angle = math.radians(38)

        max_allowable_angular_velocity = 1.0 #radians/second
        max_allowable_angle = math.radians(55)  
        dt = .01
        max_angle = 0
        max_angular_velocity = 0

        all_angles = []
        all_velocities = []

        end_effector_connections = self.end_effector_arm_connections
        l2 = self.l2

        total_time = self.max_lifetime

        if self.use_constant_period:
            organism.set_constant_period(self.constant_period)
            total_time = self.constant_period

        for i in range(int(total_time / dt)):
            time = i * dt

            results = organism.solve_sin_functions(time)

            x = results.get("x") * 0.001
            y = results.get("y") * 0.001
            z = results.get("z") * 0.001
            roll = results.get("roll")
            pitch = results.get("pitch")
            yaw = results.get("yaw")

            q = quaternion_from_euler(math.radians(roll), math.radians(pitch), math.radians(yaw))

            try:
                angles = self.solve_inverse_kinematics(translation=[x,y,z],
                                                      quaternion=q,
                                                      l2=l2,
                                                      end_effector_connections=end_effector_connections)
            except:
                self.get_logger().info(f'Failed to solve IK: {organism}')  
                return False
            
            #see if there is a new max angle
            for angle in angles:
                if abs(angle) > max_angle:
                    max_angle = abs(angle)

            for odd_arm in [angles[0], angles[2], angles[4]]:
                if odd_arm < (-1 * ground_angle):
                    hit_ground = True

            for even_arm in [angles[1], angles[3], angles[5]]:
                if even_arm > ground_angle:
                    hit_ground = True

            
            #set initial positions 
            if time == 0:
                last_angles = angles

            #calculate velocities and check for max velocity if time > 0
            else:
                velocities = [(angles[i] - last_angles[i]) / dt for i in range(len(angles))]

                for velocity in velocities:
                    if abs(velocity) > max_angular_velocity:
                        max_angular_velocity = abs(velocity)
                
                all_angles.append(angles)
                all_velocities.append(velocities)

                last_angles = angles
        

        # check if the robot even hits the ground:
        if self.check_for_hit_ground:
            if not hit_ground:

                self.get_logger().info(f'Parameters not validated due to hitting ground: {organism}')  

                return False
        
        # Check for too large angle on any arm
        
        if max_angle > max_allowable_angle:

            self.get_logger().info(f'Parameters not validated due to large angle: {organism}')  

            return False
        

        # make organism go max speed if we're using "vary_speed" setting
        if self.use_varied_speed:
            organism.vary_speed(max_angular_velocity / max_allowable_angular_velocity)   
            return True

        # check for too large velocity on any arm (only if we haven't adjusted the speed above)
        elif max_angular_velocity > max_allowable_angular_velocity:
                
            self.get_logger().info(f'Parameters not validated due to too high angular velocity \n Parameters: {organism} \n Angular velocity: {max_angular_velocity}')  
            return False
        
        # if not using varied speed and max velocity is less than allowable max
        else:
            return True
        
    def ramp_up_angles(self,angles):

        time = self.sim_time - self.sim_time_offset

        # at the begining of the organisms simulation scale up its arm angles from 0 to what it should be
        if time <= self.ramp_up_time:
            angles = [time/self.ramp_up_time * angle for angle in angles]
            return angles
        
        # at the end of each organisms simulation, scale down its arm angles to 0
        if time >= self.max_lifetime - self.ramp_up_time:
            angles = [(self.max_lifetime - time) / self.ramp_up_time * angle for angle in angles]
            return angles
        
        return angles

    def get_sim_time(self, msg):
  
        if self.running == False:
            self.sim_time_offset = self.sim_time
        self.sim_time = msg.clock.sec + (msg.clock.nanosec * 0.000000001) 
        if self.sim_time - self.sim_time_offset > self.max_lifetime:
            self.sim_time_offset = self.sim_time
            self.end_of_life()

        lifetime = self.sim_time - self.sim_time_offset
        lifetime_msg = Clock() 
        lifetime_msg.clock.sec = int(lifetime)  
        lifetime_msg.clock.nanosec = int((lifetime - lifetime_msg.clock.sec) * 1e9)
        self.lifetime_pub.publish(lifetime_msg)
    
    def get_robot_pose_from_sim(self,msg):

        if msg.header.frame_id == 'bot_world':

            if self.at_starting_position:
                self.organism_starting_position = [msg.pose.position.x,
                                                    msg.pose.position.y,
                                                    msg.pose.position.z]
                self.at_starting_position = False

            self.organism_position = [msg.pose.position.x,
                                      msg.pose.position.y,
                                      msg.pose.position.z]
            
            x = msg.pose.position.x - self.organism_starting_position[0]
            y = msg.pose.position.y - self.organism_starting_position[1]

            distance_traveled_msg = Float64()
            distance_traveled_msg.data = math.sqrt(x**2 + y**2)
            self.distance_traveled_pub.publish(distance_traveled_msg)


            sim_time = msg.header.stamp.sec + (msg.header.stamp.nanosec * 0.000000001) 
            time_alive = sim_time - self.sim_time_offset

            bot_pose_msg = String()

            data = [
                sim_time,
                time_alive,
                self.current_generation_number, 
                self.current_organism_number, 
                msg.pose.position.x,
                msg.pose.position.y,
                msg.pose.position.z,
                msg.pose.orientation.x,
                msg.pose.orientation.y,
                msg.pose.orientation.z,
                msg.pose.orientation.w,
                self.angles_from_sim[0],
                self.angles_from_sim[1],
                self.angles_from_sim[2],
                self.angles_from_sim[3],
                self.angles_from_sim[4],
                self.angles_from_sim[5],
                self.end_effector_pose.position.x,
                self.end_effector_pose.position.y,
                self.end_effector_pose.position.z,
                self.end_effector_pose.orientation.x,
                self.end_effector_pose.orientation.y,
                self.end_effector_pose.orientation.z,
                self.end_effector_pose.orientation.w
            ]

            def format_float(value):
                if isinstance(value, float):
                    return f"{value:.5f}"
                return str(value)

            formatted_data = [format_float(x) for x in data]
            bot_pose_msg.data = ','.join(formatted_data)

            self.bot_pose_pub.publish(bot_pose_msg)

    def set_end_effector_pose(self):                     #finds the pose for the end effector given the sim time and the current organisms parameters
        if self.running:
            time = self.sim_time #- self.sim_time_offset 
            current_org = self.population.organisms[self.current_organism_number]
            values_dict = current_org.solve_sin_functions(time)


            q = quaternion_from_euler(math.radians(values_dict.get('roll')),
                                      math.radians(values_dict.get('pitch')),
                                      math.radians(values_dict.get('yaw')))
            

            x = values_dict.get('x') * 0.001               #converts to meters
            y = values_dict.get('y') * 0.001
            z = values_dict.get('z') * 0.001           


            #this is used for the bot_pose_pub that saves the info for blender animations
            # this doesn't work probably because theres a delay between the goal pose and results from the sim
            self.end_effector_pose = Pose()
            self.end_effector_pose.position.x = x
            self.end_effector_pose.position.y = y
            self.end_effector_pose.position.z = z
            self.end_effector_pose.orientation.x = q[0]
            self.end_effector_pose.orientation.y = q[1]
            self.end_effector_pose.orientation.z = q[2]
            self.end_effector_pose.orientation.w = q[3]


            try:
                angles = self.solve_inverse_kinematics(translation=[x,y,z], quaternion= q)
            except:
                self.get_logger().info(f'Failed to solve IK on:  {self.population.organisms[self.current_organism_number]} \n for sim time:  {time}')  
            
            angles = self.ramp_up_angles(angles)

            try:
                msg1 = Float64()
                msg1.data = angles[0]
                
                msg2 = Float64()
                msg2.data = angles[1]

                msg3 = Float64()
                msg3.data = angles[2]

                msg4 = Float64()
                msg4.data = angles[3]

                msg5 = Float64()
                msg5.data = angles[4]

                msg6 = Float64()
                msg6.data = angles[5]

                self.pub1.publish(msg1)
                self.pub2.publish(msg2)
                self.pub3.publish(msg3)
                self.pub4.publish(msg4)
                self.pub5.publish(msg5)
                self.pub6.publish(msg6)

            except:
                self.get_logger().info(f'Failed to publish angle messages')  

            if self.publish_joint_goal_states:
              try:
                  msg = JointState()
                  msg.position = angles
                  self.joint_goal_pub.publish(msg)
              except: 
                  self.get_logger().info(f'Failed to publish angle goals JS messages')
    
    def end_of_life(self): 
        self.running = False
        distance_traveled = math.sqrt((self.organism_position[0] - self.organism_starting_position[0])**2 +
                                      (self.organism_position[1] - self.organism_starting_position[1])**2 +
                                      (self.organism_position[2] - self.organism_starting_position[2])**2)
        
        self.population.organisms[self.current_organism_number].fitness = distance_traveled

        self.at_starting_position = True

        self.get_logger().info(f'Distances traveled: {[f"{organism.fitness:.4f}"for organism in self.population.organisms]}') 

        if self.current_organism_number == (len(self.population.organisms) - 1):       #check if last organism in population
            self.new_generation()
            self.current_organism_number = 0
            self.current_generation_number += 1

            if self.use_ik_genes:
                self.end_effector_arm_connections = self.find_servo_or_connection_point_coords(self.population[self.current_organism_number][6][0] * 0.001,
                                                                                               self.population[self.current_organism_number][6][1] * 0.001,
                                                                                               end_connection_z_offset)
                self.l2 = self.population[self.current_organism_number][6][2] * 0.001

            self.get_logger().info(f'Org #{self.current_organism_number} phenotype:\n{self.population.organisms[self.current_organism_number].get_phenotype_dict()}\n')
            self.publish_pop_info()

            self.running = True
            

        else:
            self.current_organism_number += 1

            if self.use_ik_genes:
                self.end_effector_arm_connections = self.find_servo_or_connection_point_coords(self.population[self.current_organism_number][6][0] * 0.001,
                                                                                               self.population[self.current_organism_number][6][1] * 0.001,
                                                                                               end_connection_z_offset)
                self.l2 = self.population[self.current_organism_number][6][2] * 0.001

            self.get_logger().info(f'Org #{self.current_organism_number} phenotype:\n{self.population.organisms[self.current_organism_number].get_phenotype_dict()}\n')
            self.publish_pop_info()
            self.running = True
                  
    def new_generation(self):
        last_gen_number = self.population.number
        next_generation = Population(number=last_gen_number+1)

        for i, org in enumerate(self.population.organisms):
            self.get_logger().info(f'current population, org {i}: {org.genotype}')

        #selection
        selected_parents = selection(self.population, self.use_fitness_scaling)

        for i, org in enumerate(selected_parents):
            self.get_logger().info(f'selected parents, org {i}: {org.genotype}')


        #crossover
        for i in range(len(selected_parents)//2):

            parent1 = copy.deepcopy(selected_parents[i*2])
            parent2 = copy.deepcopy(selected_parents[i*2 +1])

            rand = random.random()
            if rand <= self.crossover_frequency:
                offspring = self.find_valid_offspring(parent1, parent2)

                next_generation.organisms.extend(offspring)
            else:
                next_generation.organisms.extend([parent1, parent2])
           

        for i, org in enumerate(next_generation.organisms):
            self.get_logger().info(f'offspring, org {i}: {org.genotype}')


        #mutation
        for organism in next_generation.organisms:
            self.find_valid_mutation(organism)

        for i, org in enumerate(next_generation.organisms):
            self.get_logger().info(f'mutations, org {i}: {org.genotype}')

        self.population = next_generation

        self.get_logger().info(f'Mutation completed\n New Population: {self.population}')

        return
    
    def find_valid_mutation(self, organism: Organism, attempts = None):

        if attempts == None:
            attempts = 100
        
        #validate a mutation of the organism
        for i in range(attempts):
            #use a copy of the organism, so the mutation isn't permanent
            organism_copy = copy.deepcopy(organism)
            organism_copy.mutate(mutation_frequency=self.mutation_frequency)

            if self.check_parameters(organism_copy):

                organism.genotype = organism_copy.genotype
                organism.update_phenotype()

                return
        
        self.get_logger().warn(f"Unable to produce valid mutation after {attempts} attempts")
        return

    def find_valid_offspring(self, parent1: Organism, parent2: Organism) -> List[Organism]:

        possible_cut_points = list(range(len(parent1.genotype)))
        random.shuffle(possible_cut_points)

        for cut_point in possible_cut_points:
            offspring1, offspring2 = crossover(parent1, parent2, cut_point)
            if self.check_parameters(offspring1) and self.check_parameters(offspring2):
                return [offspring1, offspring2]
        

        # Log a message if no valid offspring were found and return the parents
        self.get_logger().warn(f"Unable to produce valid offspring from given parents, returning the parents instead.")
        return [parent1, parent2]

    def solve_inverse_kinematics(self, translation, quaternion, l2 = None, end_effector_connections = None):   

        if l2 == None:
            l2 = self.l2

        if end_effector_connections == None:
            end_effector_connections = self.end_effector_arm_connections

        angles = []

        for arm_number in range(6):

            #rotate the end connection point

            rotation_matrix = R.from_quat(quaternion)
            end_connection_point = rotation_matrix.apply(end_effector_connections[arm_number])

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
            
    def find_servo_or_connection_point_coords(self, offset_1, offset_2, z_height):

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

def quaternion_from_euler(ai, aj, ak):
    ai /= 2.0
    aj /= 2.0
    ak /= 2.0
    ci = math.cos(ai)
    si = math.sin(ai)
    cj = math.cos(aj)
    sj = math.sin(aj)
    ck = math.cos(ak)
    sk = math.sin(ak)
    cc = ci*ck
    cs = ci*sk
    sc = si*ck
    ss = si*sk

    q = np.empty((4, ))
    q[0] = cj*sc - sj*cs
    q[1] = cj*ss + sj*cc
    q[2] = cj*cs - sj*sc
    q[3] = cj*cc + sj*ss

    return q
            
def main(args=None):
    rclpy.init()                        
    node = SimulationManager()      

    try:
        rclpy.spin(node)              #Starts the node and keeps it running untill it's interrupted by the user
    except KeyboardInterrupt:
        print("Terminating Node... ")
        node.destroy_node()


if __name__ == '__main__':
    main()