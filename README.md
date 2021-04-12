# Repo for control subsytem of team SeDriCa (UMIC)

Pre-requisites: 
> ROS installed in the system
> download the k.bag file using this [link](https://drive.google.com/drive/u/1/folders/1CddKUEha1fm97YP0JH-WBm0rT8SY2SPO).

### Installation:

1. Open the terminal (using ctrl+alt+T)
2. git clone https://github.com/Dikshuy/SeDriCa-Controls.git
3. cd ~/SeDriCa-Controls/controls_ws/
4. catkin_make
5. source devel/setup.bash

### Usage:

Run the following commands in different terminals:
* roscore
* rosrun Astar_ros Astar_node
* rosrun velocity_plan velocity_plan_node
* rosbag play k.bag
* rosrun Astar_ros control_node.py