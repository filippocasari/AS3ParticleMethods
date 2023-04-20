import itertools
import matplotlib.pyplot as plt
import numpy as np
from Particle import Particle as p
import math
import random
L = 30
sigma = 1
epsilon = 1
m = 1
rc = 2.5
dt = 0.01
N = int(input("Insert N: "))
#N=2

def Utruncated(U_r, U_rc, r, rc):
    return U_r-U_rc if (r <= rc) else 0


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
    #velocities = [random.uniform(-10, 10) for _ in range(N)]

    #total_momentum = sum(velocities)
    #avg_velocity = total_momentum / N
    #velocities = [v - avg_velocity for v in velocities]
    # print(velocities)
    velocities =np.zeros(N)
    #print(f"Total momentum: {total_momentum}")
    print(f"Average velocity: {sum(velocities)}")
    return velocities


velocities = return_velocities()
'''
rand_x, rand_y = [11, 11]
v_x = 0
v_y = 0
x_cell = (rand_x // rc)
y_cell = (rand_y // rc)
particle = p(rand_x, rand_y, v_x, v_y, m, x_cell, y_cell, L)
particles.append(particle)
rand_x, rand_y = [12, 12]
v_x = 0
v_y = 0
x_cell = (rand_x // rc)
y_cell = (rand_y // rc)
particle = p(rand_x, rand_y, v_x, v_y, m, x_cell, y_cell, L)
particles.append(particle)

'''
for i in range(N):
    rand_x, rand_y = random.choice(positions)
    positions.remove([rand_x, rand_y])
    v_x = (velocities[i]*math.cos(math.pi/4))
    v_y = (velocities[i]*math.sin(math.pi/4))
    x_cell = int(rand_x / rc)
    y_cell = int(rand_y / rc)
    # print(v_x, v_y)

    particle = p(rand_x, rand_y, v_x, v_y, m, x_cell, y_cell, L)
    particles.append(particle)



scatter_particles = ax.scatter([part.x for part in particles], [
    part.y for part in particles], c='red')


plt.draw()
plt.pause(1)
counter = 0

list_mom = []
potential_energy_list = []
kinetic_energy_list = []

for iter in range(10000):
    
    for index, part in enumerate(particles):
        particles[index].update(dt, rc)
    
    counter += 1
    potential_energy = 0.
    for index, part in enumerate(particles):
        particles[index].set_forces_zero()
        #print(f"forces: {particles[index].F_x}, { particles[index].F_y}")
        F_x, F_y = 0, 0
        # already_eaten = False
        i = part.i
        j = part.j
        x, y = part.x, part.y

        for ii, jj in itertools.product(range(i-1, i+1), range(j-1, j+1)):
            nearby_parts = [r for r in particles if r.i == ii %
                            N_X and r.j == jj % N_Y and r != part]
            # print(f"Nearby rabbits: {nearby_rabbits}")
            
            
            for part2 in nearby_parts:
                
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
                    
                    tg = d_y/d_x
                    angle = math.atan(tg)

                    F_x += ((48.*(d_x/r))/r**2)*(1./(r**12) - 0.5/(r**6))
                    F_y += ((48.*(d_y/r))/r**2)*(1./(r**12) - 0.5/(r**6))
                    
                    potential_energy += 4 * \
                        (1./r**12 - 1./r**6) - 4*(1./rc**12 - 1./rc**6)
                    
        particles[index].update_vel_acc(dt, F_x, F_y)
        print(f"(after) forces (particle {index}): {particles[index].F_x}, { particles[index].F_y}")
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
