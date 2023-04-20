import itertools
import matplotlib.pyplot as plt
import numpy as np
from Particle_simple import Particle as p
import math
import random
import time
L = 30
sigma = 1
epsilon = 1
m = 1
rc = 2.5
dt = 0.01
N = int(input("Insert N: "))
#N=2




N_X = int(L/rc)
N_Y = int(L/rc)
x_init = np.arange(0, L, 1)
y_init = np.arange(0, L, 1)
positions = []
for i in x_init:
    positions.extend([i, j] for j in y_init)
print(x_init)
particles = []

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
ax.grid(True, which='both', linestyle='-', linewidth=1)
ax.set_xticks(np.arange(0, L+rc, rc))
ax.set_yticks(np.arange(0, L+rc, rc))
ax.set_xlim([0, L])
ax.set_ylim([0, L])

momentum = 0
T = 0.1

def return_momentum_x(particles):
    vel_x = [part.vx for part in particles]
    return sum(vel_x)


def return_momentum_y(particles):
    v_y = [part.vy for part in particles]
    return sum(v_y)


def return_velocities():
    
    #temp = 0.1
    # return math.sqrt(3*temp/m)
    velocities = [random.uniform(-1, 1) for _ in range(N)]

    total_momentum = sum(velocities)
    avg_velocity = total_momentum / N
    velocities = [v - avg_velocity for v in velocities]
    # print(velocities)
    #velocities =np.zeros(N)
    #print(f"Total momentum: {total_momentum}")
    print(f"Average velocity: {sum(velocities)}")
    return velocities


velocities = return_velocities()

for i in range(N):
    rand_x, rand_y = random.choice(positions)
    positions.remove([rand_x, rand_y])
    v_x = np.random.randint(low=-1, high=1)*(velocities[i]*math.cos(math.pi/4))
    v_y = np.random.randint(low=-1, high=1)*(velocities[i]*math.sin(math.pi/4))
    

    particle = p(rand_x, rand_y, v_x, v_y, m, L, 0, 0)
    particles.append(particle)



scatter_particles = ax.scatter([part.x for part in particles], [
    part.y for part in particles], c='red')


plt.draw()
plt.pause(1)
counter = 0

list_mom = []
potential_energy_list = []
kinetic_energy_list = []
start = time.time()
for iter in range(200):
    
    for index, part in enumerate(particles):
        particles[index].update(dt, rc)
    
    counter += 1
    potential_energy = 0.
    for index, part in enumerate(particles):
        particles[index].set_forces_zero()
        #print(f"forces: {particles[index].F_x}, { particles[index].F_y}")
        F_x, F_y = 0, 0
        # already_eaten = False
        
        x, y = part.x, part.y

        
        for part2 in particles:
            if(part2 == part):
                continue
            d_x = x-part2.x
            if (d_x > L/2):
                d_x = d_x-L
            elif (d_x <= -L/2):
                d_x = d_x+L
            d_y = y-part2.y
            if (d_y > L/2):
                d_y = d_y-L
            elif (d_y <= -L/2):
                d_y = d_y+L
            #if (d_x == 0 and d_y == 0):
            #    d_x = np.finfo(float).eps
            #    d_y = np.finfo(float).eps

            r = math.sqrt(d_x**2 + d_y**2)
            
            if (r <= rc and r>0):

                F_x += ((48.*(d_x/r))/r**2)*(1./(r**12) - 0.5/(r**6))
                F_y += ((48.*(d_y/r))/r**2)*(1./(r**12) - 0.5/(r**6))
                
                potential_energy += 4 * \
                    (1./r**12 - 1./r**6) - 4*(1./rc**12 - 1./rc**6)
                    
        particles[index].update_vel_acc(dt, F_x, F_y)
        #print(f"(after) forces (particle {index}): {particles[index].F_x}, { particles[index].F_y}")
    potential_energy_list.append(potential_energy)
    kinetic_energy = sum([part.k_en for part in particles])
    print("delta energy: ", kinetic_energy-potential_energy)
    kinetic_energy_list.append(kinetic_energy)
    # assert(potential_energy + kinetic_energy == 0)
    velocities_x =[]
    velocities_y = []
    for index, part in enumerate(particles):
        #particles[index].set_avg_vel(return_momentum_x(
        #    particles)/(N*m), return_momentum_y(particles)/(N*m))
        #particles[index].subtract_avg_vel()
        velocities_y.append(part.vy)
        velocities_x.append(part.vx)
    
    sum_momentum = math.sqrt((sum(velocities_x)**2)+(sum(velocities_y)**2))
    list_mom.append(sum_momentum)

    
    ax.clear()
    ax.grid(True, which='both', linestyle='-', linewidth=1)
    ax.set_xticks(np.arange(0, L+rc, rc))
    ax.set_yticks(np.arange(0, L+rc, rc))
    ax.set_xlim([0, L])
    ax.set_ylim([0, L])
    scatter_particles = ax.scatter([part.x for part in particles], [
        part.y for part in particles], c='red')
    ax2.plot(list(range(iter+1)), list_mom)
    ax2.grid(True)
    ax2.set_ylabel("Momentum")
    ax2.set_xlabel("Iter")
    ax3.plot(list(range(iter+1)), kinetic_energy_list,
             c='red', label="Kinetic")
    ax3.set_ylabel("Energy")
    ax3.set_xlabel("Iter")
    ax3.plot(list(range(iter+1)), potential_energy_list,
             c='blue', label="Potential")
    ax3.plot(list(range(iter+1)), np.array(potential_energy_list)+np.array(kinetic_energy_list),\
             c='green', label="Potential+Kinetic")
    # ax3.legend().remove()
    legend = ax3.legend()
    legend.remove()

    fig.canvas.draw()
    fig.canvas.flush_events

    plt.pause(0.00001)
print(f"time execution {time.time()-start}")