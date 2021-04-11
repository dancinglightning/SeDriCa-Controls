#include<iostream>
#include<stdlib.h>
#include<vector>
#include<time.h>
#include<math.h>

#include "ros/ros.h"
#include "std_msgs/Header.h"
#include "geometry_msgs/PoseStamped.h"
#include "nav_msgs/OccupancyGrid.h"
#include "nav_msgs/Path.h"

#include <opencv2/core/core.hpp>     // opencv libraries

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;
using namespace std;

struct my_Point
{
  int x, y;
};

float width=5; //width of car in pixels,  1.5m / 0.3


ros::Publisher path_pub, pub_blown_local_map;




// float calc_gcost(float old_g,int parent_x,int parent_y,int x, int y)
// {
//   float g;
//   if(parent_x-x==1 || parent_x-x==-1)
//     { 
//       if(parent_y-y==1 || parent_y-y==-1)
//         g=old_g +14;

//       if(parent_y-y==0)
//         g=old_g +10;
//     }

//     if(parent_x-x==0)
//       g=old_g +10;    

//   return g;
// }

// float calc_gcost(node* n_pt, node* parent_n_pt)
// {
//   float g;
//   if((n_pt->parent_x)-(n_pt->x)==1 || (n_pt->parent_x)-(n_pt->x)==-1)
//     { 
//       if((n_pt->parent_y)-(n_pt->y)==1 || (n_pt->parent_y)-(n_pt->y)==-1)
//         g=parent_n_pt->g +14;

//       if((n_pt->parent_y)-(n_pt->y)==0)
//         g=parent_n_pt->g +10;
//     }

//     if((n_pt->parent_x)-(n_pt->x)==0)
//       g=parent_n_pt->g +10;    

//   return g;
// }

int index_closed=-1;   // storing the position of last filled element
int index_open=-1;

int start_x,start_y, end_x = 130 ,end_y = 70 ;


bool is_in_list(int list[][2],int i,int j,int index)
{
	for(int a=0;a<=index;a++)
	{	
		if(list[a][0]==i && list[a][1]==j)
		{
			return true;
		}
	}	
	return false;
}


class node 
{
	public:
	int parent_x,parent_y,x,y; 
	float f,g,h;
	bool is_obs,is_path;
	
	node()
	{
		h=-1;g=-1;f=-1;parent_x=-1;parent_y=-1;x=-1;y=-1;is_obs=false,is_path=false;
	}

	void calc_gcost();


	void calc_h()
	{
		h=(abs(end_x-x)+abs(end_y-y));
	}

	void calc_f()
	{
		calc_h();
    calc_gcost();

    // cout<<"\nggggggggggggg : "<<g<<"  "<<x<<"  "<<y;

		f=g+h;
		
	}

};


vector< vector<node> > n;   ///////////// DECLARING A 2D VECTOR OF NODE OBJECTS


void node::calc_gcost()
{
  if(parent_x-x==1 || parent_x-x==-1)
  { 
    if(parent_y-y==1 || parent_y-y==-1)
      g=n[parent_x+1][parent_y+1].g +14;

    if(parent_y-y==0)
      g=n[parent_x+1][parent_y+1].g +10;
  }

  if(parent_x-x==0)
    g=n[parent_x+1][parent_y+1].g +10;    
}


void SetGoal(const geometry_msgs::PoseStamped::ConstPtr& end)
{
  end_x=int((end->pose.position.x)/0.3);
  end_y=int((end->pose.position.y)/0.3);
}

float width_assign_left(int k, int a)
  { 
    float width_left=0;

   for(int j=0;j<=width;j++)
		{ if  (n[a-j+1][k+1].is_obs==false)
          { width_left= j;}

          if (n[a-j+1][k+1].is_obs==true)
           {break;}
        }

 		return width_left;
     }
       
float width_assign_right(int k, int a)
  { 
    float width_right=0;

    for(int j=0;j<=width;j++)
		{ if  (n[a+j+1][k+1].is_obs==false)
          {width_right= j;}

          if (n[a+j+1][k+1].is_obs==true)
            {break;}
        }

        return width_right;
     }


void OccupancyGridCallback(const nav_msgs::OccupancyGrid::ConstPtr& msg)

{
  if(!msg->data.empty())
  {
    int i,j,MapW,MapH;
    float mapRes;
    MapW = msg->info.width;   //  THIS MAP WIDTH AND MAP HEIGHT ARE IN PIXELS
    // cout<<msg->header.frame_id;
    // cout << MapW;
    // cout << "\n";
    MapH = msg->info.height;
    mapRes=msg->info.resolution;
    // cout<<"mapRessssss  "<<mapRes<<"\n ";
    // cout << MapH;
    // cout << "\n";
    start_x = MapW/2;
    // start_x=0;
    start_y = 0;

    // cout<<"\n start_x :"<<start_x;
    // cout<<"\n start_y :"<<start_y<<"\n";
    
    vector<my_Point> Path_nodes;
    nav_msgs::Path path;
    path.header.frame_id = "/map";


    nav_msgs::OccupancyGrid DilatedMap;
    DilatedMap.header.stamp = ros::Time::now();
    DilatedMap.header.frame_id = "/map";
    DilatedMap.info.resolution = 0.3;
    DilatedMap.info.origin.position.x = 0.0;
    DilatedMap.info.origin.position.y = 0.0;
    DilatedMap.info.origin.position.z = 0.0;
    DilatedMap.info.origin.orientation.x = 0.0;
    DilatedMap.info.origin.orientation.y = 0.0;
    DilatedMap.info.origin.orientation.z = 0.0;
    DilatedMap.info.origin.orientation.w = 1.0;
    DilatedMap.info.width = MapW;
    DilatedMap.info.height = MapH;
    DilatedMap.info.map_load_time = ros::Time::now();

    for(int i=0;i<MapH;i++)
    {
      for(int j=0;j<MapW;j++)
      {
        DilatedMap.data.push_back(0);
      }
    }

    // cout<<"tttttttttttt = "<<total_nodes;

    // node n[MapW+2][MapH+2];

    ///// CLEARING THE PREVIOUSLY STORED MAP

    n.clear();
    index_closed=-1;   //initializiing again // storing the position of last filled element
    index_open=-1;
    Path_nodes.clear();   //refreshing the path nodes ...just as path
    path.poses.clear();



    /////////// INITIALISING THE VECTOR 
    for(int i=0;i<MapH+2;i++)
    {
      vector<node> v;
      for(int j=0;j<MapW+2;j++)
      {
        node n0;
        v.push_back(n0);
      }
      n.push_back(v);
    }
	
   //  void calc_gcost()
   // {
   //  if(parent_x-x==1 || parent_x-x==-1)
   //  { 
   //    if(parent_y-y==1 || parent_y-y==-1)
   //      g=n[parent_x+1][parent_y+1].g +14;

   //    if(parent_y-y==0)
   //      g=n[parent_x+1][parent_y+1].g +10;
   //  }

   //  if(parent_x-x==0)
   //    g=n[parent_x+1][parent_y+1].g +10;    

   // }



  cv::Mat LocalMap= Mat::zeros(MapH,MapW,CV_8UC1);


   for (int i=0; i<MapH; i++)
   {
     for (int j=0; j<MapW; j++)
     { int i1=MapH-i-1;
       LocalMap.at<uchar>(i1,j) = msg->data[i*MapW+j];
     }
   }

  // cv::imshow("mat_object", LocalMap);
  // cv::waitKey(1);

  Mat dilated_mat;
  dilate(LocalMap, dilated_mat, getStructuringElement(MORPH_RECT, Size(9, 9)));

  // cv:imshow("dilated_mat", dilated_mat);
  // cv::waitKey(1);

  for (int i=0; i<MapH; i++)
  {
    for (int j=0; j<MapW; j++)
    { 
      int i1=MapH-i-1;
      DilatedMap.data[i1*MapW+j]=dilated_mat.at<uchar>(i,j);

    }
  }

    // for (int i=0; i<MapH; i++)
    // {
    //   for (int j=0; j<MapW; j++)
    //   { 
    //     // int i1=MapH-i-1;
    //     DilatedMap.data[i*MapW+j]=msg->data[i*MapW+j];
    //     cout<<DilatedMap.data[i*MapW+j];
    //     cout<<"\n\n";
    //   }
    // }



    for(i=0;i<MapH;i++)//putting values
    {
      for(int j1=MapH-1,j=0;j1>=0,j<MapW;j1--,j++)
      {
        // cout<<"\nyoooooooooooo_ "<<i<<", "<<j;
        if (DilatedMap.data[i*MapW+j]!=0)
        {
         n[j+1][i+1].is_obs=true; 
        }
      }
    }


    // cout<<"\n outtttttttt ";


    int total_nodes;
    total_nodes=MapH*MapW;

	  int closed_list[total_nodes][2];
	  int open_list[total_nodes][2];
	

	  for(int a=0;a<total_nodes;a++)                   // initializing the array
	  {	
   	  for(int b=0;b<2;b++)
   	  {	
   	  	closed_list[a][b]=-1;
   	  	open_list[a][b]=-1;
   	  }	
   	}  

    
    
    // int i,j;


   	// now making the borders with obstacles
 
   	for(i=0,j=0;j<MapH+2;j++)    
   	{	n[i][j].is_obs=true;
   	}	

   	for(i=1,j=MapH+1;i<=MapW;i++)
   	{	n[i][j].is_obs=true;
   	}	

   	for(i=MapW+1,j=0;j<MapH+2;j++)    
   	{	n[i][j].is_obs=true;
   	}	


   	for(i=1,j=0;i<=MapW;i++)
   	{	
   		n[i][j].is_obs=true;
   	}	



   	for(i=0;i<MapW+2;i++)                  //alloting coordinates to the nodes
   	{	
   		for(j=0;j<MapH+2;j++)
   		{
   			n[i][j].x=i-1;
   			n[i][j].y=j-1;
   		}	
   	}	



	  int num_obs,xo,yo;
		
	  n[start_x+1][start_y+1].parent_x=-1;
	  n[start_x+1][start_y+1].parent_y=-1;
	  n[start_x+1][start_y+1].g=0.0;
	  closed_list[0][0]=start_x;
	  closed_list[0][1]=start_y;
	  index_closed++;
	

    i=start_x;
    j=start_y;

	  int min_f,min_index;  // for finding the min f in the open list and noting its index

    while(true)
    {
    	for(int a=i-1;a<=i+1;a++)
    	{	
    		for(int b=j-1;b<=j+1;b++)
    		{	
    			if(a!=i || b!=j)
    			{	
    				  if(n[a+1][b+1].is_obs==false && is_in_list(closed_list,a,b,index_closed)==false && is_in_list(open_list, a,b,index_open)==false)
        	    {	
   		   		
        			 open_list[index_open+1][0]=a;
        			 open_list[index_open+1][1]=b;
        			 index_open++;

        			 n[a+1][b+1].parent_x=i;
        			 n[a+1][b+1].parent_y=j;

               ///////// calulating g_cost

              // if(n[a+1][b+1].parent_x-n[a+1][b+1].x==1 || n[a+1][b+1].parent_x-n[a+1][b+1].x==-1)
              // { 
              //   if(n[a+1][b+1].parent_y-n[a+1][b+1].y==1 || n[a+1][b+1].parent_y-n[a+1][b+1].y==-1)
              //     n[a+1][b+1].g=n[n[a+1][b+1].parent_x+1][n[a+1][b+1].parent_y+1].g +14;

              //   if(n[a+1][b+1].parent_y-n[a+1][b+1].y==0)
              //     n[a+1][b+1].g=n[n[a+1][b+1].parent_x+1][n[a+1][b+1].parent_y+1].g +10;
              // }

              // if(n[a+1][b+1].parent_x-n[a+1][b+1].x==0)
              //   n[a+1][b+1].g=n[n[a+1][b+1].parent_x+1][n[a+1][b+1].parent_y+1].g +10;    

              // if(i-a==1 || i-a==-1)
              // { 
              //   if(j-b==1 || j-b==-1)
              //     n[a+1][b+1].g=n[i+1][j+1].g +14;

              //   if(j-b==0)
              //     n[a+1][b+1].g=n[i+1][j+1].g +10;
              // }

              // if(i-a==0)
              //   n[a+1][b+1].g=n[i+1][j+1].g +10;   

        			 n[a+1][b+1].calc_f();

               // cout<<"\n\n ffffffffff"<<n[a+1][b+1].f; 
        
        	    }		
        	}    
        }    
    	}
        
    
        //////////////  now checking the node with the lowest  f cost in the open list....


        min_f=99999;

        for(int a=0;a<=index_open;a++)
        {
        	if(n[open_list[a][0]+1][open_list[a][1]+1].f<min_f)    //open_list[a][0] will give the x coordinate of the node
        	{
        		min_f=n[open_list[a][0]+1][open_list[a][1]+1].f;
        		min_index=a;
        	}	
        }

    
        ////moving to the node with min f
        i=open_list[min_index][0];
        j=open_list[min_index][1];

        ////adding the node with minimun f to the closed list

        closed_list[index_closed+1][0]=i;// or open_list[min_index][0];   // for x 
        closed_list[index_closed+1][1]=j;//or open_list[min_index][1];	  // for y
        index_closed++;

        if(i==end_x && j==end_y)
        	break;

        ////now removing that node from the open list

        for(int a=min_index;a<=index_open-1;a++)
        {
        	open_list[a][0]=open_list[a+1][0];
        	open_list[a][1]=open_list[a+1][1];
        }
        index_open-=1;       //one element has decreased in open list	
    }


    /////we have reached the destination...
    /////so now tracing back...
   
    // publisher(n, path_pub, poses)
    int a,b,temp;
    a=end_x;b=end_y;

    //cout<<"\n\n Paaaaaaaaaaaaath \n\n";

    while(a>=0 && b>=0)
     { 
      struct my_Point p={a,b};
      //cout<<"\n ("<<p.x<<","<<p.y<<")";
      Path_nodes.push_back(p);
      n[a+1][b+1].is_path=true;
      temp=a; 
      a=n[a+1][b+1].parent_x;
      b=n[temp+1][b+1].parent_y; 

     }
    // for(int i=0;i<path.poses.size();i++)
      // cout<<Path_nodes[i].x<<Path_nodes[i].y<<"\n";

    reverse(Path_nodes.begin(),Path_nodes.end());



      
    for(int i=0; i<Path_nodes.size(); i++)
      {
        // cout<<"there is dataaaaaaa  "<<Path_nodes.size();
        geometry_msgs::PoseStamped vertex;
        vertex.header.stamp=ros::Time::now();
        vertex.header.frame_id="/map";
        vertex.pose.position.x= double((Path_nodes[i].x)*0.3);  //  BECAUSE MAP RESOLUTION IS 0.3
        vertex.pose.position.y= double((Path_nodes[i].y)*0.3);
        vertex.pose.position.z= 0.0;

        vertex.pose.orientation.x= float(width_assign_left(Path_nodes[i].y, Path_nodes[i].x))*0.3;
        vertex.pose.orientation.y= float(width_assign_right(Path_nodes[i].y, Path_nodes[i].x))*0.3;
        vertex.pose.orientation.z= 0.0;
        vertex.pose.orientation.w= 0.0;

        // cout<<"\n"<<vertex.pose;
      
        path.poses.push_back(vertex);
      }
          

    //    ////////////////////
    // cout<<"\n\n\n\n\n\n NEWWWWWWWWWWW MAPPPPPPPPPP \n\n\n\n";
    // for(int i=MapH-1;i>=0;i--)
    // {
    //   for(int j=0;j<133;j++)
    //     {if(n[j+1][i+1].is_path)
    //       cout<<"1";
    //     else if(n[j+1][i+1].is_obs)
    //       cout<<"2";
    //     else if(!n[j+1][i+1].is_obs)
    //       cout<<"0";}
    //   cout<<"\n";
    // }
    ////////////////////

    pub_blown_local_map.publish(DilatedMap);
    path_pub.publish(path);

  }

  // else
  //   cout<<"\n msg is emptyyyyyyyy";

	
}


int main(int argc, char **argv)
{
  ros::init(argc,argv,"Astar_node");
  ros::NodeHandle n;
  ros::Subscriber goal_sub = n.subscribe<geometry_msgs::PoseStamped>("/move_base_simple/goal", 1, SetGoal);

  path_pub= n.advertise< nav_msgs::Path >("/A_star_path", 1); 
  ros::Subscriber map_sub = n.subscribe<nav_msgs::OccupancyGrid>("/scan/local_map", 1, OccupancyGridCallback); 
  //ros::Subscriber map_sub = n.subscribe<nav_msgs::OccupancyGrid>("/occ_map", 1, OccupancyGridCallback);

  pub_blown_local_map = n.advertise<nav_msgs::OccupancyGrid>("/dilated_map",1);


  /*while(ros::ok())
  {
    ros::spinOnce();
    ros::Rate rate(10);
    rate.sleep();
  }  

  return 0;*/

  ros::Rate RosLoopRate(20);     // 15 in RRT*
  while(ros::ok())
  {
    
    ros::spinOnce();//check for incoming messages
    RosLoopRate.sleep(); 

  }
  return 0; 
}