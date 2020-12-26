import numpy as np  

class Bicycle:
	#dynamic model
	#front-wheel drive 
	#slip angle = steering angle
	# u = [a, ω] control inputs 
	# ξ = [x, y, v, θ, ϕ, δ]

	def __init__(self):
        self.v_x = 0
        self.v_y = 0
        self.theta = 0  #orientation 
        # self.phi = 0    #orientation speed
        self.delta = 0  #steering angle
        self.density = 1.2 #densoty of air
        self.C_drag = 0.28  #drag coefficient 
        self.length =       #front length of car
        self.breadth =       #front breadth of car
        self.A_front = 0.85 * self.length * self.breadth

        self.J = 0 #rotational momentum 

        self.lf = 1 #distance b/w front wheels and CG
        self.lr = 1 #distance b/w rear wheels and CG
        self.m = 200 #mass of vehicle
        self.Cy = 0.1 #lateral tire stiffness
        self.g = 9.81  #gravitation constant
        self.fr = 0.01 #rolling resistance coefficient
        self.beta_x = 0.04 #road tilt in x
        self.beta_y =     # raod tilt in y

        self.sample_time = 0.01 


    def update(self, a, w):
        F_drag = 0.5 * self.density * self.C_drag * self.A_front * self.v_x**2 * np.sign(self.v_x)   #drag force
        F_rolling = self.fr * self.m * self.g * np.min(1, v_x) * np.sign(v_x)                        #rolling force
        F_slope = self.m * self.g * np.sin(self.beta_x) * np.sign(v_x)                               #sloping force
        F_tilt = self.m * self.g * np.sin(self.beta_y)                                               #tilt force
        #following equations need to be improvised
        m*a_x = Fx1*cos(delta)-Fy1*sin(delta)+Fx2*cos(delta)-Fy2*sin(delta)+Fx3+Fx4-F_drag-F_rolling-F_slope
        m*a_y = Fx1*sin(delta) - Fy1*cos(delta) + Fx2*sin(delta) - Fy2*cos(delta) - F_tilt
        I_z * alpha_z =  (Fx1*sin(delta) + Fy1*cos(delta))*lf + (Fx2*sin(delta) + Fy2*cos(delta))*lf + (Fy3 + Fy4)*lr
     	# Fyf = self.Cy * (self.delta - (self.La*self.phi)/self.v)   #front tire lateral force
     	# Fyr = self.Cy * (self.Lb*self.phi)/self.v                  #rear tire lateral force

     	#differential equations
   #   	xc_dot = self.v * np.cos(self.theta)
   #   	yc_dot = self.v * np.sin(self.theta)
   #   	v_dot = a * np.cos(self.delta) - (2/self.m) * self.Fyf * np.sin(self.delta)
   #      theta_dot = self.phi
   #      phi_dot = (1/self.J) * (self.La*(self.m * a * np.sin(self.delta) + 2* Fyf*np.cos(self.delta)) - 2*self.Lb*Fyr)
 		# delta_dot = w 
        a_x = v_x_dot + v_y*omega_z
        a_y = v_y_dot + v_x*omega_z
        alpha_z = omega_z_dot
        omega_z = theta_z_dot

        #balance over roll axis 
        (Ix + m*d_roll**2)*alpha_x + m*d_roll*alpha_y + (K_rollf + K_rollr - m*g*d_roll)*theta_roll + (D_rollf + D_rollr)*omega_x = 0
        omega_x = theta_roll_dot
        alpha_x = omega_x_dot

        #balance over pitch axis
        (Iy + m*d_pitch**2)*alpha_y - m*d_pitch*a_x + (K_pitch + m*g*d_pitch)*theta_pitch + D*pitch*omega_y = 0
        omega_y = theta_pitch_dot
        alpha_y = omega_y_dot

        #vertical motion 
        z = COG_z + d_pitch*(cos(theta_pitch)-1) + d_roll*(cos(theta_roll)-1)
        v_z = z_dot
        a_z = v_z_dot

        #vel vector vs longitudnal axis
        theta_vf = (y_dot + lf*omega_z) / v_x
        theta_vr = (y_dot - lr*omega_z) / v_x

        #geometric vel
        v_x1 = v_x - omega_z * df/2
        v_x1 = v_x + omega_z * df/2
        v_x3 = v_x - omega_z * dr/2
        v_x4 = v_x + omega_z * dr/2



 		#update equations (not completed)
        self.v += v_dot * self.sample_time
        self.theta += theta_dot * self.sample_time
        # self.phi += phi_dot * self.sample_time
        self.delta += delta_dot * self.sample_time

        pass
