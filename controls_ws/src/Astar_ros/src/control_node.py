#!/usr/bin/env python
# -*- coding: utf-8 -*-
import rospy
from std_msgs.msg import Float64MultiArray
from std_msgs.msg import Float32
# from geometry_msgs.msg import PoseStamped
# from nav_msgs.msg import OccupancyGrid
from nav_msgs.msg import Path
import numpy as np
import sys
import os
import do_mpc
from do_mpc.data import save_results, load_results
from casadi import *
import matplotlib.pyplot as plt
# import matplotlib as mpl
import pickle
sys.path.append('../../')


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


def derivative_bezier(a,b,c,d):
    return lambda t: np.power(1-t,2)*a*(-3)+3*b*(np.power(1-t,2)-2*t*(1-t))+3*c*(2*t*(1-t)-np.power(t,2))+3*d*np.power(t,2)


def derivative_list(points):
    A,B=get_bezier_coef(points)
    return [derivative_bezier(points[i],A[i],B[i],points[i+1]) for i in range(len(points)-1)]


# evalute each cubic curve on the range [0, 1] sliced in n points
def evaluate_bezier(points, n):
    curves = get_bezier_cubic(points)
    return np.array([fun(t) for fun in curves for t in np.linspace(0, 1, n)])


def trajectory_gen(points):
    A,B=get_bezier_coef(points)
    return [get_cubic(points[i],A[i],B[i],points[i+1]) for i in range(len(points)-1)]


def control(x_0,i,x,y,vel,acc,steer_rate,points,curves,derivatives,velocities,acc_pub,brake_pub,steer_pub):
    model_type='discrete'
    model=do_mpc.model.Model(model_type)
    J=1000
    La=1
    Lb=1
    m=200
    Cy=0.1
    t_s=0.01  # sample time
    N=70
    k=0.1
    fn=curves[i]
    d=derivatives[i]
    vmax_i = 1.5
    # vmax_i=max(velocities[i],velocities[i+1])
    # state variables
    psi=model.set_variable(var_type='_x',var_name='psi',shape=(1,1))
    xc=model.set_variable(var_type='_x',var_name='xc',shape=(1,1))
    yc=model.set_variable(var_type='_x',var_name='yc',shape=(1,1))
    v=model.set_variable(var_type='_x',var_name='v',shape=(1,1))
    theta=model.set_variable(var_type='_x',var_name='theta',shape=(1,1))
    phi=model.set_variable(var_type='_x',var_name='phi',shape=(1,1))
    delta=model.set_variable(var_type='_x',var_name='delta',shape=(1,1))
    a_s=model.set_variable(var_type='_x',var_name='a_s',shape=(1,1))
    w_s=model.set_variable(var_type='_x',var_name='w_s',shape=(1,1))
    # k=model.set_variable(var_type='_p',var_name='k',shape=(1,1))
    Fyf=Cy*(delta-(La*phi)/v)
    Fyr=(Cy*Lb*phi)/v
    # control inputs
    a=model.set_variable(var_type='_u',var_name='a',shape=(1,1))
    omega=model.set_variable(var_type='_u',var_name='omega',shape=(1,1))
    model.set_expression(expr_name='cost', expr=sum1((xc-fn(psi)[0])**2+100*(yc-fn(psi)[1])**2+100*theta**2+200*(np.tan(theta)-d(psi)[1]/d(psi)[0])**2 + a_s**2+w_s**2+(v-0.9*vmax_i)**2))
    # model.set_expression(expr_name='cost', expr=sum1((xc-fn(psi)[0])**2+100*(yc-fn(psi)[1])**2+100*theta**2+200*(np.tan(theta)-d(psi)[1]/d(psi)[0])**2
                                                    # +a_s**2+w_s**2)+(v-0.8*velocities[i+1])**2)
    state_now=vertcat(psi,xc, yc, v, theta, phi, delta,a_s,w_s)
    # B=t_s*vertcat((0.9*vmax_i)/((d(psi)[0])**2+(d(psi)[1])**2)**0.5,v*np.cos(theta), v*np.sin(theta), a* np.cos(delta)-(2.0/m)*Fyf*np.sin(delta), phi,
    #             (1.0/J)*(La*(m*a*np.sin(delta)+2*Fyf*np.cos(delta))-2*Lb*Fyr), omega,(1/t_s)*(a-a_s),(1/t_s)*(omega-w_s))
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
    mpc.set_rterm(a=0.1)
    mpc.set_rterm(omega=0.1)
    # 'tube' from path planning
    # mpc.bounds['lower','_x','xc']=x_i[0]-1e-19
    mpc.bounds['lower','_x','yc']=min(points[i+1][1],x_0[2])-1.5
    mpc.bounds['lower','_x','v']=0  # max reverse speed in m/s
    mpc.bounds['lower','_x','theta']=-np.pi/2-1e-2
    # mpc.bounds['upper','_x','xc']=target+5
    mpc.bounds['upper','_x','yc']=max(points[i+1][1],x_0[2])+1.5
    mpc.bounds['upper','_x','theta']=np.pi/2 + 1e-2
    mpc.bounds['upper','_x','v']=vmax_i
    # mpc.bounds['upper','_x','phi']=50
    '''
    mpc.bounds['upper','_x','delta']=50
    mpc.bounds['lower','_u','a']=-10
    mpc.bounds['lower','_u','omega']=-10
    mpc.bounds['upper','_u','a']=10
    mpc.bounds['upper','_u','omega']=10
    mpc.bounds['lower','_x','a_s']=-10
    mpc.bounds['lower','_x','w_s']=-10
    mpc.bounds['upper','_x','a_s']=10
    mpc.bounds['upper','_x','w_s']=10
    '''
    mpc.bounds['lower','_u','omega']=-5
    mpc.bounds['upper','_u','omega']=5
    mpc.bounds['lower','_u','a']=-1
    mpc.bounds['upper','_u','a']=1
    mpc.setup()
    # estimator=do_mpc.estimator.StateFeedback(model)
    simulator = do_mpc.simulator.Simulator(model)
    simulator.set_param(t_step=t_s)
    '''p_template=simulator.get_p_template()
    def p_fun(t_now):
        p_template['k']=0.1
        return p_template
    simulator.set_p_fun(p_fun)'''
    simulator.setup()
    mpc.x0 = x_0
    simulator.x0 = x_0
    mpc.u0['a']=x_0[-2]
    mpc.u0['omega']=x_0[-1]
    mpc.set_initial_guess()
    mpc.reset_history()
    N_u=N

    # Customizing Matplotlib:
    # mpl.rcParams['font.size'] = 18
    # mpl.rcParams['lines.linewidth'] = 3
    # mpl.rcParams['axes.grid'] = True

    # mpc_graphics = do_mpc.graphics.Graphics(mpc.data)
    # sim_graphics = do_mpc.graphics.Graphics(simulator.data)

    # # We just want to create the plot and not show it right now. This "inline magic" supresses the output.
    # fig, ax = plt.subplots(2, sharex=True, figsize=(16,9))
    # fig.align_ylabels()

    # for g in [sim_graphics, mpc_graphics]:
    # # Plot the angle positions (phi_1, phi_2, phi_2) on the first axis:
    # g.add_line(var_type='_x', var_name='xc', axis=ax[0])
    # # g.add_line(var_type='_x', var_name='phi_2', axis=ax[0])
    # # g.add_line(var_type='_x', var_name='phi_3', axis=ax[0])

    # # Plot the set motor positions (phi_m_1_set, phi_m_2_set) on the second axis:
    # g.add_line(var_type='_x', var_name='yc', axis=ax[1])
    # # g.add_line(var_type='_u', var_name='phi_m_2_set', axis=ax[1])

    # ax[0].set_ylabel('X')
    # ax[1].set_ylabel('Y')
    # ax[1].set_xlabel('time [s]')

    # sim_graphics.plot_results()
    # # Reset the limits on all axes in graphic to show the data.
    # sim_graphics.reset_axes()
    # # Show the figure:
    # import time
    # plt.savefig(str(time.time())+'.png')

    simulator.reset_history()
    for _ in range(N_u):
        u0=mpc.make_step(x_0)
        if u0[0][0]>=0:
            acc_pub.publish(u0[0][0])
        else:
            brake_pub.publish((-1)*u0[0][0])
        steer_pub.publish(u0[1][0]*t_s)
        x_0=simulator.make_step(u0)
        # with open('control_outputs.csv',mode='a') as op_file:
        #     op=csv.writer(op_file,delimiter=',')
        #     op.writerow([u0[0][0],u0[1][0]])

        if x_0[1]>=points[i+1][0]:
            break
        # x_0=estimator.make_step(y_0)

    save_results([mpc, simulator])
    with open('results/results.pkl', 'rb') as f:
        results = pickle.load(f)

    print('')
    print('')
    print('='*100)
    print(results['mpc']['_x'])
    print('='*100)
    print('')
    print('')

    os.remove("results/results.pkl")

    acc=vertcat(acc,simulator.data['_u','a',-1])
    x=vertcat(x,simulator.data['_x','xc',-1])
    y=vertcat(y,simulator.data['_x','yc',-1])
    vel=vertcat(vel,simulator.data['_x','v',-1])
    steer_rate=vertcat(steer_rate,simulator.data['_u','omega',-1])
    # z=vertcat(z,simulator.data['_x','theta',-1])

    if x_0[1]>=points[i+1][0]:
        return x_0,x,y,vel,acc,steer_rate
    else:
        # x_0=simulator.data['_x'][-1]
        return control(x_0,i,x,y,vel,acc,steer_rate,points,curves,derivatives,velocities,acc_pub,brake_pub,steer_pub)


velocities=[1.5]
path_points=[]
c=1
packed=[velocities,c]


def v_callback(v_arr):
    print("here")
#     velocities.append(v_arr.data[0])
    # '''l=len(v_arr.data)
    #             packed[0]=packed[0][:packed[1]]
    #             for i in range(l):
    #                 packed[0].append(v_arr.data[i])
    #             packed[1]+=1'''


x_initial=[[0],[0],[0],[1.5],[0],[0],[0],[0],[0]]


def path_callback(path):
    vel=packed[0]
    path_repeat=False
    x_0=np.zeros((len(x_initial),1))
    for i in range(len(x_initial)):
        x_0[i][0]=x_initial[i][0]
    l = len(path.poses)
    now_points=np.zeros((l,2))
    if len(path_points):
        if path_points[0][0]==path.poses[0].pose.position.y:
            print("###repetition detected###")
            path_repeat=True
            path_sub.unregister()
            print("#######################################path unregistered#######################################")
    if not path_repeat:
        for j in range(l):
            now_points[j][0]=path.poses[j].pose.position.y-path.poses[0].pose.position.y
            now_points[j][1]=path.poses[j].pose.position.x-path.poses[0].pose.position.x
            path_points.append([path.poses[j].pose.position.y,path.poses[j].pose.position.x])
        x_=np.array([x_0[1]])
        y_=np.array([x_0[2]])
        v=np.array([x_0[3]])
        a=np.array([x_0[7]])
        s_r=np.array([x_0[-1]])
        first,rest=now_points[:50,:],now_points[49:,:]
        while first.shape[0]==50:
            bcurves=trajectory_gen(first)
            derivatives=derivative_list(first)
            for i in range(len(first)-1):
                (x_0,x_,y_,v,a,s_r)=control(x_0,i,x_,y_,v,a,s_r,first,bcurves,derivatives,vel,acc_pub,brake_pub,steer_pub)
                x_0[0]=0
            if rest.shape[0]>50:
                first,rest=rest[:50,:],rest[49:,:]
            else:
                break
        bcurves=trajectory_gen(rest)
        derivatives=derivative_list(rest)
        for i in range(len(rest)-1):
            (x_0,x_,y_,v,a,s_r)=control(x_0,i,x_,y_,v,a,s_r,rest,bcurves,derivatives,vel,acc_pub,brake_pub,steer_pub)
            x_0[0]=0
        for i in range(len(x_initial)):
            x_initial[i][0]=x_0[i][0]
        fig,ax=plt.subplots(2,2)
        ax[0][0].plot(x_,y_)
        ax[0][0].scatter(now_points[:,0],now_points[:,1])
        ax[0][1].plot(x_,v)
        ax[0][1].plot(now_points[:,0],vel[:l],'ro')
        # opt=velocities[:l]
        for i in range(len(opt)):
            opt[i]=0.8*opt[i]
        ax[0][1].plot(now_points[:,0],opt,'bo')
        ax[1][0].plot(x_,a)
        ax[1][1].plot(x_,s_r)
        ax[0][0].set(xlabel='x',ylabel='y')
        ax[0][1].set(xlabel='x',ylabel='v')
        ax[1][0].set(xlabel='x',ylabel='a')
        ax[1][1].set(xlabel='x',ylabel='steer rate')
        plt.show()


rospy.init_node('control_node', anonymous=True)
acc_pub = rospy.Publisher('acceleration', Float32, queue_size=10)
brake_pub = rospy.Publisher('brake', Float32, queue_size=10)
steer_pub = rospy.Publisher('steer', Float32, queue_size=10)
rate=rospy.Rate(10)
vel_sub=rospy.Subscriber("/velocity_plan",Float64MultiArray,v_callback)
path_sub=rospy.Subscriber("/A_star_path",Path,path_callback)
rospy.spin()
