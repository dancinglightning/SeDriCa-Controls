import numpy as np
import sys
sys.path.append('../../')
import do_mpc
from casadi import *
import matplotlib.pyplot as plt
N_ref=300
x_ref=np.linspace(0,20,N_ref+1)
y_ref=0.15*np.sin(x_ref/5)
model_type='discrete'
J=500
La=1
Lb=1
m=200
Cy=0.1
t_s=0.01 #sample time
N=30
def controller_fn(model_i,x_0,target_x,y_upper,y_lower,i):
    mpc=do_mpc.controller.MPC(model_i)
    setup_mpc = {
        'n_horizon': N,
        't_step': t_s,
        'n_robust': 0,
        'state_discretization':'discrete',
        'store_full_solution': True,
    }
    mpc.set_param(**setup_mpc)
    mterm = model_i.aux['cost']
    lterm = model_i.aux['cost']
    mpc.set_objective(mterm=mterm, lterm=lterm)
    mpc.set_rterm(a=0.1)
    mpc.set_rterm(w=0.1)
    # 'tube' from path planning
    mpc.bounds['lower','_x','xc']=x_0[0]-1e-4
    mpc.bounds['lower','_x','yc']=y_lower
    mpc.bounds['lower','_x','v']=0 #max reverse speed in m/s
    mpc.bounds['lower','_x','theta']=-50
    mpc.bounds['lower','_x','phi']=-50
    mpc.bounds['lower','_x','delta']=-50
    mpc.bounds['upper','_x','xc']=target_x+0.1
    mpc.bounds['upper','_x','yc']=y_upper
    mpc.bounds['upper','_x','v']=20 #max forward speed in m/s
    mpc.bounds['upper','_x','theta']=50
    mpc.bounds['upper','_x','phi']=50
    mpc.bounds['upper','_x','delta']=50
    mpc.bounds['lower','_u','a']=-10
    mpc.bounds['lower','_u','w']=-10
    mpc.bounds['upper','_u','a']=10
    mpc.bounds['upper','_u','w']=10
    mpc.bounds['lower','_x','a_s']=-10
    mpc.bounds['lower','_x','w_s']=-10
    mpc.bounds['upper','_x','a_s']=10
    mpc.bounds['upper','_x','w_s']=10
    mpc.setup()
    mpc.x0 = x_0
    if i>0:
        mpc.u0['a']=x_0[-2]
        mpc.u0['w']=x_0[-1]
    mpc.set_initial_guess()
    mpc.reset_history()
    u0=mpc.make_step(x_0)
    return u0
def update(i,x,y,x_ref,y_ref,x_0):
    model=do_mpc.model.Model(model_type)
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
    w=model.set_variable(var_type='_u',var_name='w',shape=(1,1))
    state_now=vertcat(xc, yc, v, theta, phi, delta,a_s,w_s)
    B=t_s*vertcat(v*np.cos(theta), v*np.sin(theta), a_s* np.cos(delta)-(2.0/m)*Fyf*np.sin(delta), phi,
                    (1.0/J)*(La*(m*a*np.sin(delta)+2*Fyf*np.cos(delta))-2*Lb*Fyr), w,(1/t_s)*(a-a_s),(1/t_s)*(w-w_s))
    state_next=state_now + B
    model.set_rhs('xc',state_next[0])
    model.set_rhs('yc',state_next[1])
    model.set_rhs('v',state_next[2])
    model.set_rhs('theta',state_next[3])
    model.set_rhs('phi',state_next[4])
    model.set_rhs('delta',state_next[5])
    model.set_rhs('a_s',state_next[6])
    model.set_rhs('w_s',state_next[7])
    target_x=x_ref[i+1]
    target_y=y_ref[i+1]
    if x_0[0]>=target_x:
        target_x=x_0[0]+x_ref[i+1]-x_ref[i]
        #target_y=0.15*np.sin(target_x/5)
    ##doubt
    x_range=np.linspace(x_0[0],target_x,100)
    y_range=0.15*np.sin(x_range/5)
    (y_upper,y_lower)=(max(y_range)+1e-1,min(y_range)-1e-1)
    ##doubt
    model.set_expression(expr_name='cost', expr=sum1((xc-target_x)**2+(yc-target_y)**2)
                        +(a_s)**2+w_s**2)
    model.setup()
    simulator = do_mpc.simulator.Simulator(model)
    simulator.set_param(t_step = t_s)
    simulator.x0 = x_0
    simulator.setup()
    for j in range(N):
        u0_i=controller_fn(model,x_0,target_x,y_upper,y_lower,i)
        x_0=simulator.make_step(u0_i)
        if x_0[0]>=target_x:
            breakS
    x=vertcat(x,simulator.data['_x','xc',-1])
    y=vertcat(y,simulator.data['_x','yc',-1])
    return x,y,x_0
x_0=np.array([[0],[0],[0.001],[(np.pi/4)*(0.03)],[0],[0],[0],[0]])
#initial condition on theta
x=np.zeros((1,1))
y=np.zeros((1,1))
for i in range(N_ref):
    x,y,x_0=update(i,x,y,x_ref,y_ref,x_0)
plt.plot(x,y)
plt.plot(x_ref,y_ref-1e-1)
plt.plot(x_ref,y_ref+1e-1)