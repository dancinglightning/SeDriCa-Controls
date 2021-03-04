#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include<time.h>
#include<math.h>
//#include<Python.h>

#include "ros/ros.h"
#include "std_msgs/Header.h"
#include "std_msgs/Float32.h"
#include "std_msgs/Float32MultiArray.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/OccupancyGrid.h"
#include "nav_msgs/Path.h"

using namespace std;

float x_c=0;
float y_c=0;
vector<float> x_coordinates;
vector<float> y_coordinates;
std_msgs::Float32MultiArray msg_x;
std_msgs::Float32MultiArray msg_y;
/*void Callback(const nav_msgs::Path::ConstPtr& msg)
{
  
  ros::Publisher x_pub = n.advertise<std_msgs::Float32>("x_c", 1000);
  ros::Publisher y_pub = n.advertise<std_msgs::Float32>("y_c", 1000);
  if(!msg->poses.empty()){
    
    x_coordinates.clear();
    y_coordinates.clear();
    int data = msg -> poses.size();
    for(int i=0; i<data; i++){
      y_c = msg->poses[i].pose.position.x;
      x_c = msg->poses[i].pose.position.y;
      std_msgs::Float32 msg_x;
      std_msgs::Float32 msg_y;
      msg_x.data=x_c;
      msg_y.data=y_c;
      x_pub.publish(msg_x);
      y_pub.publish(msg_y);
     //fstream fout;
     //fout.open("coordinates.csv",ios::out | ios::app);
     //fout<<x_c<<","<<y_c<<"\n";
//      ROS_INFO("[%f]", x_c);
//      ROS_INFO("[%f]", y_c);
      x_coordinates.push_back(x_c);
      y_coordinates.push_back(y_c);
      std::cout<<x_c<<' '<<y_c<<'\n';
    }
    //system("python3 mpc.py");
  }
}*/

class SubscribeAndPublish
{
public:
  SubscribeAndPublish()
  {
    //Topic you want to publish
    x_pub = n.advertise<std_msgs::Float32MultiArray>("x_c_vector", 1);
    y_pub = n.advertise<std_msgs::Float32MultiArray>("y_c_vector", 1);

    //Topic you want to subscribe
    sub = n.subscribe("/A_star_path", 1, &SubscribeAndPublish::callback, this);
  }

  void callback(const nav_msgs::Path::ConstPtr& msg)
  {
    if(!msg->poses.empty()){
    
    x_coordinates.clear();
    y_coordinates.clear();
    msg_x.data.clear();
    msg_y.data.clear();
    int data = msg -> poses.size();

    for(int i=0; i<data; i++){
      y_c = msg->poses[i].pose.position.x;
      x_c = msg->poses[i].pose.position.y;
      
      msg_x.data.push_back(x_c);
      msg_y.data.push_back(y_c);
      if (msg_x.data.size()==20 and msg_y.data.size()==20){
         x_pub.publish(msg_x);
         y_pub.publish(msg_y);
         msg_x.data.clear();
         msg_y.data.clear();}
     //fstream fout;
     //fout.open("coordinates.csv",ios::out | ios::app);
     //fout<<x_c<<","<<y_c<<"\n";
//      ROS_INFO("[%f]", x_c);
//      ROS_INFO("[%f]", y_c);
      x_coordinates.push_back(x_c);
      y_coordinates.push_back(y_c);
      std::cout<<x_c<<' '<<y_c<<'\n';
    }
   
    //system("python3 mpc.py");
  }
  }

private:
  ros::NodeHandle n; 
  ros::Publisher x_pub;
  ros::Publisher y_pub;
  ros::Subscriber sub;

};//End of class SubscribeAndPublish

int main(int argc, char **argv)
{
  //Initiate ROS
  ros::init(argc, argv, "coordinates");

  //Create an object of class SubscribeAndPublish that will take care of everything
  SubscribeAndPublish SAPObject;

  ros::spin();

  return 0;
}


/*int main(int argc, char *argv[])
{
  ros::init(argc, argv, "coordinates");

  
 //for(int i=0; i < x_coordinates.size(); i++)
   //{myFile<< x_coordinates.at(i) << ','<<y_coordinates.at(i)<<'\n';}
//std::cout<<x_c<<' '<<y_c<<std::endl;
  //runPython(argc, argv);
  ros::spin();

  return 0;
}*/
