#launch file for seeing the urdf model move around with the gui

from launch import LaunchDescription
from launch_ros.actions import Node
from launch.actions import ExecuteProcess
from launch.actions import IncludeLaunchDescription
from launch_ros.substitutions import FindPackageShare
from launch.substitutions import PathJoinSubstitution
from launch.launch_description_sources import PythonLaunchDescriptionSource
import xacro

def generate_launch_description():
    return LaunchDescription([      

    #launches the simulation in gazebo
        ExecuteProcess(
            cmd=['ign', 'gazebo', PathJoinSubstitution([FindPackageShare("genetic_algorithm"), "sdf","my_bot.sdf"]), '-r'],   #'-r', '-s' 
            output="screen"
        ),

    #Node for managing the simulation
        Node(
            package= "genetic_algorithm",
            executable= "ign_gazebo_sim.py",
            name= "ign_gz_sim_Node"
        ),
    
    #node for bridging topics between ros and gazebo
        Node(
            package= "ros_gz_bridge",
            executable= "parameter_bridge",
            name= "ros_gz_bridge_node",
            parameters=[{'config_file': PathJoinSubstitution([FindPackageShare("genetic_algorithm"), "config", "ros_gz_bridge.yaml"])}]
        )  

    ])
