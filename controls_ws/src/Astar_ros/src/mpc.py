import numpy as np
import do_mpc
import sys
sys.path.append('../../')
from casadi import *

model_type = 'continuous'  # either 'discrete' or 'continuous'
model = do_mpc.model.Model(model_type)

# model variables

# states
x = model.set_variable(var_type='_x',var_name='x',shape=(1,1))  # x-position of the vehicle in the reference frame
y = model.set_variable(var_type='_x',var_name='y',shape=(1,1))  # y-position of the vehicle in the reference frame
theta = model.set_variable(var_type='_x',var_name='theta',shape=(1,1))  # orientation 
delta = model.set_variable(var_type='_x',var_name='delta',shape=(1,1))  # steering angle
x_dot = model.set_variable(var_type="_x", var_name='x_dot', shape=(1,1))
y_dot = model.set_variable(var_type="_x", var_name='y_dot', shape=(1,1))
theta_dot = model.set_variable(var_type="_x", var_name='theta_dot', shape=(1,1))

#inputs
omega=model.set_variable(var_type='_u',var_name='omega',shape=(1,1))  # theta_dot
v = model.set_variable(var_type='_u',var_name='v',shape=(1,1))  # velocity 

# model parameters
L = model.set_variable('parameter', 'L')  # length of the car

L = 1

'''
equations:

x_dot_e = v*cos(theta_e) + ye*(v_d/L)*tan(delta_d) - v_d
y_dot_e = v*sin(theta_e) - xe*(v_d/L)*tan(delta_d)
theta_e = (v/L)*tan(delta) - (v_d/L)*tan(delta_d) 

'''

model.set_rhs('x_dot', x_dot)
model.set_rhs('y_dot', y_dot)
model.set_rhs('theta_dot', theta_dot)

x_dot = v*cos(theta_e) 
y_dot = v*sin(theta_e)
theta = (v/L)*tan(delta) 

model.setup()

mpc = do_mpc.controller.MPC(model)

setup_mpc = {
    'n_horizon': 20,
    't_step': 0.1,
    'n_robust': 1,
    'store_full_solution': True,
}

mpc.set_param(**setup_mpc)

mterm = x**2 + y**2 + theta**2
lterm = x**2 + y**2 + theta**2

mpc.set_objective(mterm=mterm, lterm=lterm)

# Lower bounds on states:
mpc.bounds['lower','_x', 'x'] = None
mpc.bounds['lower','_x', 'y'] = None
mpc.bounds['lower','_x', 'theta'] = -(np.pi)/4
# Upper bounds on states
mpc.bounds['upper','_x', 'x'] = None
mpc.bounds['upper','_x', 'y'] = None
mpc.bounds['upper','_x', 'theta'] = -(np.pi)/4

mpc.setup()

simulator = do_mpc.simulator.Simulator(model)

simulator.set_param(t_step = 0.1)

simulator.setup()

x0 = np.array([0, 0, 0, 0, 0, 0, 0]).reshape(-1,1)

simulator.x0 = x0
mpc.x0 = x0
mpc.set_initial_guess()

u0 = np.zeros((2,1))
for i in range(200):
    simulator.make_step(u0)

u0 = mpc.make_step(x0)

simulator.reset_history()
simulator.x0 = x0
mpc.reset_history()

# capture
for i in range(20):
    u0 = mpc.make_step(x0)
    x0 = simulator.make_step(u0)

