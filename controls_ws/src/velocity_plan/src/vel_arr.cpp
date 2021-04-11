#include <bits/stdc++.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include "ros/ros.h"
#include <nav_msgs/Path.h>
#include <std_msgs/Float64MultiArray.h>
#include <std_msgs/Bool.h>
#include <std_msgs/Float64.h>
using namespace std;

double Vin  = 0;                             //initial velocity , IDK why
double Vout = 0;                             //Ohh shit, here we go again

float TLightDis;                             //distance of traffic light
float TLightDir;                             //direction of traffic light
bool IsOutVel;                               //Some shit

double WtCurrRadii = 0.7;                    //weight of current radii over previous one

double Acc = 1.0;                            //Acceleration
double Decc = 1.5;                           //IDK
double MaxAcc = 1.7;                         //max centri-acceleration of car; per 1m of radii
double VelMax = 4;                           //max velocity of car
double VelZebra = 2.77;                      //max velocity of car while on zebra

int VelDCurrCan = 0;
float VelDCurr = VelDCurrCan*5/18;

std::vector<float> Zebra;                    //signal received from zebra channel, mostly distance
std::vector<bool> PathZ;                     //storing values of z-co of path

std_msgs::Float64MultiArray Vel,VelBack;     //calculated velocity backup

ros::Publisher pub_vel_arr;

void OutVelCallback(const std_msgs::Bool::ConstPtr& x)
{
	IsOutVel = x->data;
  	if(IsOutVel)
	  	Vout = 5.55;
  	else 
  		Vout = 0;
}
void TLightDisCallback(const std_msgs::Float64::ConstPtr& tra)
{
	TLightDis = (tra->data);
}
void TLightDirCallback(const std_msgs::Float64::ConstPtr& tra)
{
	TLightDir = (tra->data);
}
void ZebraCallback(const std_msgs::Float64MultiArray::ConstPtr& zeb){
	//cout<<"zebra_init";
	Zebra.clear();
	for(int i=0;i<zeb->data.size();i++)
		Zebra.push_back(zeb->data[i]);
}

vector<double> vmax(double c[], int n)
{
	vector<double> v_m;
	v_m.push_back(Vin);
	for(int i=0; i<n; i++)
	{
		if(c[i]==-1)
			v_m.push_back(0);
		else if(c[i]==-2)
			v_m.push_back(VelZebra);
		else if(c[i]>=2000)
			v_m.push_back(VelMax);
		else
			v_m.push_back(sqrt(MaxAcc*c[i]));
	}
	v_m.push_back(Vout);
	return v_m;
}
vector<vector<double> > segment(vector<double> v)
{
	int n=0,ctr = 0;
	vector<vector<double> > x;
	for(int i=0;i<v.size()-1;i++)
	{
		vector<double> row;
		while(i<v.size()-1 && v[i]>=v[i+1])
		{
			row.push_back(v[i]);
			i++;
			ctr = 1;
		}
		while(i<v.size()-1 && v[i]<=v[i+1] && ctr!=1)
		{
			row.push_back(v[i]);
			i++;
			ctr = 2;
		}
		ctr=0;
		row.push_back(v[i]);
		x.push_back(row);
	}
	return x;
}

void smooth(vector<double> &v_m, vector<double> &distance)
{
	for(int i=1;i<v_m.size();i++)
	{
		if((v_m[i-1]-v_m[i])>1)
		{
			for(int j=i; j>=0; j--)
			{
				double val2 = sqrt( v_m[i]*v_m[i]+2*Acc*(distance[i]-distance[j]));	
				v_m[j]=min(val2,v_m[j]);}
		}
	}
}

void PathCallback(const nav_msgs::Path::ConstPtr& msg){
	if(!msg->poses.empty())
	{
		PathZ.clear();
		int NumData = msg->poses.size();

		double RadiiPath[NumData-2];
		Vel.data.push_back(Vin); //initial condition

		vector<double> PathDis;
		vector<double> VelArr;
		vector<vector<double> > SegVelArr;
		for(int i=NumData-2,k=0; i>=1; i--,k++)
		{
			double x1 = msg->poses[i-1].pose.position.x ,x2 = msg->poses[i].pose.position.x, x3 = msg->poses[i+1].pose.position.x;
			double y1 = msg->poses[i-1].pose.position.y ,y2 = msg->poses[i].pose.position.y, y3 = msg->poses[i+1].pose.position.y;
			PathZ.push_back(msg->poses[i].pose.position.z);

			double a=(x1*(y2-y3)-y1*(x2-x3)+x2*y3-y2*x3);
			double b=((x1*x1+y1*y1)*(y3-y2)+(x2*x2+y2*y2)*(y1-y3)+(x3*x3+y3*y3)*(y2-y1));
			double c=((x1*x1+y1*y1)*(x2-x3)+(x2*x2+y2*y2)*(x3-x1)+(x3*x3+y3*y3)*(x1-x2));
			double d=((x1*x1+y1*y1)*(x3*y2-x2*y3) + (x2*x2+y2*y2)*(x1*y3-x3*y1) + (x3*x3+y3*y3)*(x2*y1-x1*y2));

			if(a!=0)
				RadiiPath[k]=sqrt((b*b+c*c-4*a*d)/(4*a*a));
			if(a==0)
				RadiiPath[k]=2000; //max radii is 2000 ----- note for whole script
			PathDis.push_back(0.0);
			for(int j=1; j<NumData; j++)
				PathDis.push_back(PathDis[j-1]+sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))));
		}
		for(int i=1;i<NumData-2;i++)
			RadiiPath[i] = WtCurrRadii*RadiiPath[i] + (1-WtCurrRadii)*RadiiPath[i-1];
		for(int i=0;i<NumData-2;i++)
			if(RadiiPath[i]>(VelMax*VelMax/MaxAcc))
				RadiiPath[i]=2000;

		//Traffic Light Detection and Distance
		int TLightPos = NumData-1;
		for(int i=0;i<NumData-1;i++)
		{
			if(TLightDis>=(PathDis[i]+7) && TLightDis<=(PathDis[i]+8) && TLightDir==1)
			{
				RadiiPath[i]=-1;
				TLightPos = i;
				break;
			}
		}
		if(TLightDis>=0 && TLightDis<=7)
			TLightPos = 0;
		for(int i=TLightPos;i<NumData-2;i++)
			RadiiPath[i] = -1;

		//Zebra Crossing
		for(int i=0;i<Zebra.size();i++)
		{
			if(Zebra[i]!=-1)
			{
				for(int j=0;j<NumData-2;j++)
				{
					if(PathDis[j]==Zebra[i])
						RadiiPath[j]=-2;
				}
			}
		}

		VelArr = vmax(RadiiPath,NumData-2);
		SegVelArr = segment(VelArr);
		smooth(VelArr,PathDis);

		int a=0;
		int b=0;
		for(int i=0; i<SegVelArr.size();i++)
		{
			for(int j=0; j<SegVelArr[i].size();j++)
			{
				if(!(i==0 && j==0))
				{
					double val1 = sqrt(Vel.data[b]*Vel.data[b]+2*Acc*(PathDis[a]-PathDis[b]));
					double val2 = sqrt(VelArr[b+SegVelArr[i].size()-1]*VelArr[b+SegVelArr[i].size()-1]+2*Decc*(PathDis[b+SegVelArr[i].size()-1]-PathDis[a]) ) ;
					double temp = min((double)VelArr[a],(double)val1);
					Vel.data.push_back(min((double)temp, (double)val2));
				}
				a++;
			}
			b = b + SegVelArr[i].size();
			Vel.data[b] = Vel.data[b-1];
		}
		VelDCurr = VelDCurrCan*5/18;
		Vin = VelDCurr;
		std::reverse(Vel.data.begin(),Vel.data.end());
		pub_vel_arr.publish(Vel);
		VelBack = Vel;
		Vel.data.clear();
	}
	else
	{
		for(int i=0; i<VelBack.data.size(); i++)
		{
			VelBack.data[i] = 0;
		}
		pub_vel_arr.publish(VelBack);
	}
}

int main(int argc, char **argv)
{
	//In the x vector the fields are as follows: 1. x 2. y 3. traffic sense (+1 or -1) 4.zebra sense 5.Distance to traffic lights (distance or -1)
	ros::init(argc, argv, "vel_arr");
	ros::NodeHandle RosNodeH;

	TLightDir = -1;
	TLightDis = -1;
	Zebra.push_back(-1);

	ros::Subscriber sub_output_vel  = RosNodeH.subscribe("/output_vel",1,OutVelCallback);
	ros::Subscriber sub_traffic_dis = RosNodeH.subscribe("/traffic_distance",1,TLightDisCallback);
	ros::Subscriber sub_traffic_dir = RosNodeH.subscribe("/traffic_direction",1,TLightDirCallback);
	ros::Subscriber sub_zebra       = RosNodeH.subscribe("/zebra",1,ZebraCallback);
	ros::Subscriber sub_path        = RosNodeH.subscribe("/A_star_path",1,PathCallback); 
		
	pub_vel_arr = RosNodeH.advertise<std_msgs::Float64MultiArray>("/velocity_array", 1);

	ros::Rate loop_rate(200);

	while (ros::ok())
	{
		ifstream myfile_read("../../../scripts/can_output.txt"); 
		if(myfile_read.is_open())
		{
			myfile_read>>VelDCurrCan;
			myfile_read.close();
		} 
		// for(int j = 0;j<v_copy.data.size();j++)
		// 	cout<<j<<" "<<v_copy.data[j]<<"  "<<rev[j]<<endl;
		// cout<<" ///////////////////////////////////////////////////////////// "<<endl; 
		ros::spinOnce();
		loop_rate.sleep();
	}
	return 0;
}