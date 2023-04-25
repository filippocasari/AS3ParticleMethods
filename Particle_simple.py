import math
import numpy as np

class Particle:
    def __init__(self, x, y, vx, vy, m,  L, i, j, thermostat_temp):

        self.x = x

        self.y = y
        self.thermostat_temp = thermostat_temp
        self.vx = vx
        self.vy = vy
        self.m = m
        self.L = L
        self.a_x = 0.
        self.a_y = 0.
        self.F_x = 0.
        self.F_y = 0.
        self.i = i
        self.j =j
        # self.avg_vel_x =0
        # self.avg_vel_y =0
        self.k_en = 0.
        self.compute_k_energy()
    # update position

    def update(self, dt, rc):
        self.y += (self.vy*dt)+self.a_y*dt*dt*0.5
        self.x += (self.vx*dt) + self.a_x*dt*dt*0.5
        self.x = self.x % self.L
        self.y = self.y % self.L
        self.i = int(self.x//rc)
        self.j = int(self.y//rc)
        

    def update_vel_acc(self, dt, F_x, F_y):
        
        
        self.F_x = F_x
        self.F_y = F_y
        new_acc_x = self.F_x/self.m
        new_acc_y = self.F_y/self.m
        self.vx += (self.a_x+new_acc_x)*dt*0.5
        self.vy += (self.a_y+new_acc_y)*dt*0.5
        self.a_x = new_acc_x
        self.a_y = new_acc_y
        self.compute_k_energy()

    def set_forces(self, F_x, F_y):
        self.F_x = F_x
        self.F_y = F_y


    def compute_k_energy(self):
        self.k_en = 0.5*self.m*((self.vx**2)+(self.vy**2))


    def set_forces_zero(self):
        self.F_x = 0.
        self.F_y = 0.
        #self.a_x =0.
        #self.a_y=0.
    def apply_thermostat_scaling(self, temp):
        
        gamma = np.sqrt(1+(0.0025*((self.thermostat_temp/temp)-1)))
        self.vx *= gamma
        self.vy *= gamma