#!/usr/bin/env python
# -*- coding: utf-8 -*-
import rospy
from std_msgs.msg import Float64MultiArray
from std_msgs.msg import Float32
from nav_msgs.msg import Path
import numpy as np
import sys
import os
sys.path.append('../../')
import do_mpc
from casadi import *
import matplotlib.pyplot as plt
import csv
from do_mpc.data import save_results, load_results

def get_bezier_coef(points):
    n = len(points) - 1

    # build coefficents matrix
    C = 4 * np.identity(n)
    np.fill_diagonal(C[1:], 1)
    np.fill_diagonal(C[:, 1:], 1)
    C[0, 0] = 2
    C[n - 1, n - 1] = 7
    C[n - 1, n - 2] = 2

    # build points vector
    P = [2 * (2 * points[i] + points[i + 1]) for i in range(n)]
    P[0] = points[0] + 2 * points[1]
    P[n - 1] = 8 * points[n - 1] + points[n]

    # solve system, find a & b
    A = np.linalg.solve(C, P)
    B = [0] * n
    for i in range(n - 1):
        B[i] = 2 * points[i + 1] - A[i + 1]
    B[n - 1] = (A[n - 1] + points[n]) / 2

    return A, B

# returns the general Bezier cubic formula given 4 control points
def get_cubic(a, b, c, d):
    return lambda t: np.power(1 - t, 3)*a + 3*np.power(1 - t, 2)*t*b + 3*(1 - t)*np.power(t, 2)*c + np.power(t, 3)*d

# return one cubic curve for each consecutive points
def get_bezier_cubic(points):
    A, B = get_bezier_coef(points)
    return [
        get_cubic(points[i], A[i], B[i], points[i + 1])
        for i in range(len(points) - 1)
    ]

# return derivative of the cubic function with given coefficients
def derivative_bezier(a,b,c,d):
    return lambda t: np.power(1-t,2)*a*(-3)+3*b*(np.power(1-t,2)-2*t*(1-t))+3*c*(2*t*(1-t)-np.power(t,2))+3*d*np.power(t,2)

# return derivative for each consecutive points
def derivative_list(points):
    A,B=get_bezier_coef(points)
    return [derivative_bezier(points[i],A[i],B[i],points[i+1]) for i in range(len(points)-1)]

# return list of lambda functions for each consecutive points
def trajectory_gen(points):
    A,B=get_bezier_coef(points)
    return [get_cubic(points[i],A[i],B[i],points[i+1]) for i in range(len(points)-1)]

# function to update state and publish acc, steer, brake
def control(x_0,x,y,vel,acc,steer_rate,curves,derivatives,points,velocities,acc_pub,brake_pub,steer_pub,i):
	model_type='discrete'
	model=do_mpc.model.Model(model_type)
	J=1000 #rotational momentum
	La=1
	Lb=1
	m=200
	Cy=0.1 # lateral tire stiffness
	t_s=0.01 #sample time
	N=60
	k=0.1
	fn=curves[i] #follow i-th cubic curve
	d=derivatives[i] #derivative of i-th cubic fn
	vmax_i=velocities[i] #max velocity for between i-th and i+1-th point
	#state variables
	psi=model.set_variable(var_type='_x',var_name='psi',shape=(1,1))
	xc=model.set_variable(var_type='_x',var_name='xc',shape=(1,1))
	yc=model.set_variable(var_type='_x',var_name='yc',shape=(1,1))
	v=model.set_variable(var_type='_x',var_name='v',shape=(1,1))
	theta=model.set_variable(var_type='_x',var_name='theta',shape=(1,1))
	phi=model.set_variable(var_type='_x',var_name='phi',shape=(1,1))
	delta=model.set_variable(var_type='_x',var_name='delta',shape=(1,1))
	a_s=model.set_variable(var_type='_x',var_name='a_s',shape=(1,1))
	w_s=model.set_variable(var_type='_x',var_name='w_s',shape=(1,1))
	Fyf=Cy*(delta-(La*phi)/v)
	Fyr=(Cy*Lb*phi)/v
	#control inputs
	a=model.set_variable(var_type='_u',var_name='a',shape=(1,1))
	omega=model.set_variable(var_type='_u',var_name='omega',shape=(1,1))
	model.set_expression(expr_name='cost', expr=sum1((xc-fn(psi)[0])**2+100*(yc-fn(psi)[1])**2+100*theta**2+200*(np.tan(theta)-d(psi)[1]/d(psi)[0])**2
                                                     +a_s**2+w_s**2)+(v-0.8*vmax_i)**2)
	state_now=vertcat(psi,xc, yc, v, theta, phi, delta,a_s,w_s)
	B=t_s*vertcat(k,v*np.cos(theta), v*np.sin(theta), a* np.cos(delta)-(2.0/m)*Fyf*np.sin(delta), phi,
                  (1.0/J)*(La*(m*a*np.sin(delta)+2*Fyf*np.cos(delta))-2*Lb*Fyr), omega,(1/t_s)*(a-a_s),(1/t_s)*(omega-w_s))
	state_next=state_now + B
	model.set_rhs('psi',state_next[0])
	model.set_rhs('xc',state_next[1])
	model.set_rhs('yc',state_next[2])
	model.set_rhs('v',state_next[3])
	model.set_rhs('theta',state_next[4])
	model.set_rhs('phi',state_next[5])
	model.set_rhs('delta',state_next[6])
	model.set_rhs('a_s',state_next[7])
	model.set_rhs('w_s',state_next[8])
	model.setup()
	mpc=do_mpc.controller.MPC(model)
	setup_mpc = {
	'n_horizon': N,
	't_step': t_s,
	'n_robust': 1,
	'state_discretization':'discrete',
	'store_full_solution': True,
	}
	mpc.set_param(**setup_mpc)
	mterm = model.aux['cost']
	lterm = model.aux['cost']
	mpc.set_objective(mterm=mterm, lterm=lterm)
	mpc.set_rterm(a=0.01)
	mpc.set_rterm(omega=0.01)
	mpc.bounds['lower','_x','yc']=min(points[1][1],x_0[2])-1.5
	mpc.bounds['lower','_x','v']=0 #max reverse speed in m/s
	mpc.bounds['lower','_x','theta']=-np.pi/2-1e-2
	mpc.bounds['upper','_x','yc']=max(points[1][1],x_0[2])+1.5
	mpc.bounds['upper','_x','theta']=np.pi/2 + 1e-2
	mpc.bounds['upper','_x','v']=vmax_i
	mpc.bounds['lower','_u','omega']=-5
	mpc.bounds['upper','_u','omega']=5
	mpc.bounds['lower','_u','a']=-1
	mpc.bounds['upper','_u','a']=1
	mpc.setup()
	simulator = do_mpc.simulator.Simulator(model)
	simulator.set_param(t_step = t_s)
	simulator.setup()
	mpc.x0 = x_0
	simulator.x0 = x_0
	mpc.u0['a']=x_0[-2]
	mpc.u0['omega']=x_0[-1]
	mpc.set_initial_guess()
	mpc.reset_history()
	for i in range(50):
		u0=mpc.make_step(x_0)
		if u0[0][0]>=0:
			acc_pub.publish(u0[0][0])
		else:
			brake_pub.publish((-1)*u0[0][0])
		steer_pub.publish(u0[1][0]*t_s)
		x_0=simulator.make_step(u0)
	acc=vertcat(acc,simulator.data['_u','a',-1])
	x=vertcat(x,simulator.data['_x','xc',-1])
	y=vertcat(y,simulator.data['_x','yc',-1])
	vel=vertcat(vel,simulator.data['_x','v',-1])
	steer_rate=vertcat(steer_rate,simulator.data['_u','omega',-1])
	return x_0,x,y,vel,acc,steer_rate

def v_callback(v_arr):
	for i in range(len(v_arr.data)):
		velocities.append(v_arr.data[i])
	global x_0
	x_0[0]=0
	x_0[1]=points[0][0]
	x_0[2]=points[0][1]
	x=np.array(x_0[1])
	y=np.array(x_0[2])
	v=np.array(x_0[3])
	acc=np.array(x_0[7])
	s_r=np.array(x_0[-1])
	bcurves=trajectory_gen(points[:20,:])
	derivatives=derivative_list(points[:20,:])
	fig,ax=plt.subplots(2,1)
	ax[0].set(xlabel='x',ylabel='y')
	ax[1].set(xlabel='x',ylabel='y')
	for i in range(10):
		x_0,x,y,v,acc,s_r=control(x_0,x,y,v,acc,s_r,bcurves,derivatives,points,velocities,acc_pub,brake_pub,steer_pub,i)
		x_0[0]=0
	try:
		ax[0].plot(x,y)
		ax[0].plot(points[:10,0],points[:10,1],'bo')
		ax[1].plot(x,v)
		ax[1].plot(points[:10,0],velocities[:10],'ro')
	except:
		ax[0].plot(points[:10,0],points[:10,1])
	plt.show()


'''def v_callback(v_arr,x_0):
	#with open('vel_new.csv',mode='a') as op_file:
		#op=csv.writer(op_file)
	#for i in range(len(v_arr.data)):
		#velocities.append(v_arr.data[i])
		#	op.writerow([v_arr.data[i]])
		#op.writerow(['##############'])
	velocities=[]
	with open('subscribed_vel.csv',mode='r') as ip_file:
		csv_reader=csv.reader(ip_file,delimiter=',')
		i=0
		for row in csv_reader:
			if row[0][0]=='#':
				break
			velocities.append(float(row[0]))
	x=np.array(x_0[1])
	y=np.array(x_0[2])
	v=np.array(x_0[3])
	acc=np.array(x_0[7])
	s_r=np.array(x_0[-1])
	first,rest=points[:20,:],points[19:,:] # get first 20 points ("first"), and store remaining in "rest"
	fig,ax=plt.subplots(2,1)
	ax[0].set(xlabel='x',ylabel='y')
	ax[1].set(xlabel='x',ylabel='v')
	j=1
	while(True):
		bcurves=trajectory_gen(first)
		derivatives=derivative_list(first)
		for i in range(first.shape[0]-1):
			x_0,x,y,v,acc,s_r=control(x_0,x,y,v,acc,s_r,bcurves,derivatives,first,velocities,acc_pub,brake_pub,steer_pub,i)
			x_0[0]=0
		try:
			ax[0].plot(x,y)
			ax[0].plot(points[:20*j,0],points[:20*j,1],'ro')
			ax[1].plot(x,v)
			ax[1].scatter(points[:,0],velocities)
		except:
			ax[0].plot(x,y)
		plt.show()
		fig.clear()
		if(rest.shape[0]<=20):
			break
		first,rest=rest[:20,:],rest[19:,:] #get next 20 points to "first"
		j+=1
	print(rest.shape)
	bcurves=trajectory_gen(rest)
	derivatives=derivative_list(rest)
	for i in range(rest.shape[0]-1):
		x_0,x,y,v,acc,s_r=control(x_0,x,y,v,acc,s_r,bcurves,derivatives,rest,velocities,acc_pub,brake_pub,steer_pub,i)
		x_0[0]=0
	ax[0].plot(x,y)
	ax[0].plot(points[:,0],points[:,1],'bo')
	ax[1].plot(x,v)
	ax[1].plot(points[:,0],velocities)
	plt.show()'''



def path_callback(path):
	global points
	points=np.zeros((len(path.poses),2))
	for i in range(len(path.poses)):
		points[i][0]=path.poses[i].pose.position.y-path.poses[0].pose.position.y
		points[i][1]=path.poses[i].pose.position.x-path.poses[0].pose.position.x
	'''if points is None:
					points=new_points
				else:
					if points.all()==new_points.all():
						return
					points=new_points'''
	'''with open('coordinates.csv',mode='a') as op_file:
					op=csv.writer(op_fil	e)
					for i in range(points.shape[0]):
						op.writerow([points[i][0],points[i][1]])
					op.writerow(["########"])'''

		

x_0=np.array([[0],[0],[0],[0.01],[0],[0],[0],[0],[0]])
velocities=[]
points=None
rospy.init_node('debug_node',anonymous=True)
acc_pub = rospy.Publisher('acceleration', Float32, queue_size=10)
brake_pub = rospy.Publisher('brake', Float32, queue_size=10)
steer_pub = rospy.Publisher('steer', Float32, queue_size=10)
rate=rospy.Rate(10)
vel_sub=rospy.Subscriber("/velocity_plan",Float64MultiArray,v_callback,queue_size=1)
path_sub=rospy.Subscriber("/A_star_path",Path,path_callback,queue_size=1)
rospy.spin()
