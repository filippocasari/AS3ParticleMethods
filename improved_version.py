import itertools
import matplotlib.pyplot as plt
import numpy as np
from Particle_simple import Particle as p
import math
import random
import time
import matplotlib.collections as mc
L = 30
sigma = 1
epsilon = 1
m = 1
rc = 2.5
dt = 0.01
N = int(input("Insert N: "))

N_X = int(L/rc)
N_Y = int(L/rc)

cell_list = [[[] for _ in range(N_Y)] for _ in range(N_X)]
x_init = np.arange(0, L, 1.)
y_init = np.arange(0, L, 1.)
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


def return_momentum_x(particles):
    vel_x = [part.vx for part in particles]
    return sum(vel_x)


def return_momentum_y(particles):
    v_y = [part.vy for part in particles]
    return sum(v_y)


init_temp = 0.1


def return_velocities():
    global init_temp
    init_total_kinetic_energy = 1.5*N*init_temp

    individual_kin_energy = init_total_kinetic_energy/N
    velocities = np.full(N, np.sqrt(2*individual_kin_energy/m))
    # velocities = [random.uniform(-10, 10) for _ in range(N)]
    print("Initial kinetics energy : ", init_total_kinetic_energy)
    print("Individual kinetics energy : ", individual_kin_energy)
    print(f"Average velocity: {np.mean(velocities)}")
    print(f"total system velocity {np.sum(velocities)}")
    return velocities


def compute_init_potential_energy(particles, rc):

    potential_energy = 0.
    equal_obj = 0
    for part in particles.copy():

        x, y = part.x, part.y

        for part2 in particles.copy():
            if (part != part2):
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
                # if (d_x == 0 and d_y == 0):
                #    d_x = np.finfo(float).eps
                #    d_y = np.finfo(float).eps

                r = np.sqrt(d_x**2 + d_y**2)

                if (r <= rc):

                    potential_energy += 4 * \
                        (1./r**12 - 1./r**6) - 4*((1./rc**12) - (1./rc**6))
            else:
                equal_obj += 1
    print(f"Equal objects: {equal_obj}")
    return potential_energy


velocities = return_velocities()
# positions[250+i]#
for i in range(N):
    alpha = 2 * math.pi * random.random()
    rand_x, rand_y = random.choice(positions)  
    positions.remove([rand_x, rand_y])
    v_x = velocities[i]*math.cos(alpha)
    v_y = velocities[i]*math.sin(alpha)
    x_cell = int(rand_x//rc)
    y_cell = int(rand_y//rc)

    particle = p(rand_x, rand_y, v_x, v_y, m, L, x_cell, y_cell)
    particles.append(particle)
    cell_list[x_cell][y_cell].append(particle)

points_whole_ax = 5 * 0.8 * 72
radius = 0.5
points_radius = 2 * radius / 1.0 * points_whole_ax
# particles=np.array(particles)
scatter_particles = ax.scatter([part.x for part in particles], [
    part.y for part in particles], c='red', s=0.5)


plt.draw()
plt.pause(1)


list_mom = []

temperatures = [init_temp]
start = time.time()
init_energy = 0
init_kinetic_energy = sum(part.k_en for part in particles)
print("Init kinetic energy: ", init_kinetic_energy)
init_potential_energy = compute_init_potential_energy(particles, rc)
total_energy = []
potential_energy_list = [init_potential_energy]
kinetic_energy_list = [init_kinetic_energy]
# single_part_U = []
for iter in range(300):

    for part in particles:
        cell_list[part.i][part.j].remove(part)

        part.update(dt, rc)

        cell_list[int(part.x//rc)][int(part.y//rc)].append(part)

    potential_energy = 0.
    # single_part_U_iter = 0
    for part in particles:
        part.set_forces_zero()
        F_x, F_y = 0, 0

        x, y = part.x, part.y
        cell_x, cell_y = int(x//rc), int(y//rc)
        for i, j in itertools.product(range(-1, 2), range(-1, 2)):
            for part2 in cell_list[(cell_x+i) % N_X][(cell_y+j) % N_Y]:
                if part2 != part:
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

                    r = np.sqrt(d_x**2 + d_y**2)

                    if (r <= rc and r > 0):

                        F_x += ((48.*(d_x/r))/r**2)*(1./(r**12) - 0.5/(r**6))
                        F_y += ((48.*(d_y/r))/r**2)*(1./(r**12) - 0.5/(r**6))
                        # print(f"potential energy: {potential_energy}")
                        # potential_energy += 4.*((1./(r**12)) - 1./(r**6))
                        potential_energy += 4. * \
                            ((1./(r**12)) - 1./(r**6)) - \
                            (4. * ((1./(rc**12)) - (1./(rc**6))))
                        '''if(part == particles[0]):
                            single_part_U_iter += 4. * \
                            ((1./(r**12)) - 1./(r**6))
                            
                            single_part_U.append(single_part_U_iter)'''
        part.update_vel_acc(dt, F_x, F_y)

    potential_energy_list.append(potential_energy)
    kinetic_energy = sum(part.k_en for part in particles)

    kinetic_energy_list.append(kinetic_energy)
    velocities_x = []
    velocities_y = []
    for part in particles:
        velocities_y.append(part.vy)
        velocities_x.append(part.vx)

    sum_momentum = np.sqrt((np.sum(velocities_x)**2) +
                           (np.sum(velocities_y)**2))*m
    list_mom.append(sum_momentum)

    ax.clear()
    ax.grid(True, which='both', linestyle='-', linewidth=1)
    ax.set_xticks(np.arange(0, L+rc, rc))
    ax.set_yticks(np.arange(0, L+rc, rc))
    ax.set_xlim([0, L])
    ax.set_ylim([0, L])
    # scatter_particles = ax.scatter([part.x for part in particles], [
    #    part.y for part in particles], c='red', s=30)
    sizes = [20. for _ in range(N)]
    xy = [[part.x, part.y] for part in particles]
    patches = [plt.Circle(center, size) for center, size in zip(xy, sizes)]
    collection = mc.CircleCollection(
        sizes, offsets=xy, transOffset=ax.transData, color='green')
    ax.add_collection(collection)
    ax2.plot(list(range(iter+1)), list_mom)
    ax2.grid(True)
    ax2.set_ylabel("Momentum")
    ax2.set_xlabel("Iter")

    ax3.plot(list(range(len(kinetic_energy_list))), kinetic_energy_list,
             c='red', label="Kinetic")
    ax3.set_ylabel("Energy")
    ax3.set_xlabel("Iter")
    ax3.plot(list(range(len(potential_energy_list))), potential_energy_list,
             c='blue', label="Potential")
    total_energy.append(kinetic_energy + potential_energy)
    ax3.plot(list(range(len(total_energy))), total_energy,
             c='green', label="Potential+Kinetic")
    ax3.plot(list(range(iter+1)), np.full(iter+1,
             potential_energy_list[0]+kinetic_energy_list[0]))
    # ax3.plot(list(range(iter+1)), np.full(iter+1,
    #        potential_energy_list[0]+kinetic_energy_list[0]+1/np.sqrt(N)))
    # ax3.plot(list(range(iter+1)), np.full(iter+1,
    #         potential_energy_list[0]+kinetic_energy_list[0]-1/np.sqrt(N)))
    temp = 2.*kinetic_energy/(3. * N)
    temperatures.append(temp)
    ax4.plot(list(range(len(temperatures))), temperatures,
             c='blue', label="Temp")
    ax4.set_ylabel("Temperature")
    ax4.set_xlabel("Iter")

    # ax3.legend().remove()
    legend = ax3.legend()
    legend.remove()

    fig.canvas.draw()
    fig.canvas.flush_events

    plt.pause(0.0000001)
fig.savefig(f"images/A_{N}.jpg")
print(f"time execution {time.time()-start}")
print("std temperature: ", np.std(temperatures))
print("fluctuations temperature order ", 1/np.sqrt(N))
print("std total energy: ", np.std(total_energy))
print("fluctuations total energy ", 1/np.sqrt(N))
# fig, ax = plt.subplots()
# ax.plot(list(range(len(single_part_U))), single_part_U)
# plt.grid(True)
# plt.show()
