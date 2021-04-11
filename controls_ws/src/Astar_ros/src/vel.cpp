#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include<time.h>
#include<math.h>

#include "ros/ros.h"
#include "std_msgs/Header.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/OccupancyGrid.h"
#include "nav_msgs/Path.h"
#include <std_msgs/Float64MultiArray.h>

using namespace std;

void Callback(const std_msgs::Float64MultiArray::ConstPtr& x)
{
  int len = x->data.size();
  double arr[len];
  for(int i=0; i<len; i++)
    {
      arr[i] = x->data[i];
    }
  double max = 0;
  for(int i=0; i<len; i++)
  {
    if(arr[i]>max){
      max = arr[i];
    }
  }
  ROS_INFO("[%f]", max);
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "vel_subscriber");

  ros::NodeHandle n;

  ros::Subscriber sub = n.subscribe("/velocity_array", 1, Callback);

  ros::spin();

  return 0;
}