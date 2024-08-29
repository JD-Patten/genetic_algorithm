
# genetic_algorithm
This is a ROS2 Package for setting up and running a genetic algorithm with the purpose of making a stewart platform walk. It works For stewart platforms that have arms composed of two links.

This package was made for the robot in this [youtube video](https://www.youtube.com/watch?v=Ly8sEk-lyek), but has the goal of being easily customizeable and expandable. 

## Need Help /  Want to Discuss
Join the [Discord](https://discord.gg/4GBbg6FE)!

If you have reccomendations, want to help out with this project, or have any questions, I'll be happy to discuss more in the discord!

## HOW IT WORKS
The main launch files starts 3 things: the simulation, the node to controll the sim, and a node to create a bridge for topics between the sim and the controller. 

The ignition gazebo simulation is set up with an sdf file called my_bot.sdf which defines the physics, the robot model, and the rest of the world.

The ROS2 node ran with the ign_gazebo_sim.py script is used for taking in all of the needed data from the simulation, sending joint state commands to the simulation to move the robot around, and do all the Genetic algorithm operations.

