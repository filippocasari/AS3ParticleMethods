#include <fstream>
#include <iostream>
#include <thread>
#include <algorithm>
#include <mutex>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <tuple>
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <netinet/in.h>
#include "Particle.h"
#include <string>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <czmq.h>
#define GLFW_INCLUDE_NONE
const int WIDTH = 800;
const int HEIGHT = 600;
using namespace std;

double* return_velocities(double init_temp, int N, double m){
    
        double init_total_kinetic_energy = 1.5*((double) N) *init_temp;

        const double individual_kin_energy = init_total_kinetic_energy/(double) N;
        //velocities = np.full(N, np.sqrt(2*individual_kin_energy/m));
        double* velocities = new double[N];
        for(int i=0; i< N; ++i){
            velocities[i]=sqrt(2.0*individual_kin_energy/m);
        }
        printf("Initial kinetics energy : %f\n", init_total_kinetic_energy);
        printf("Individual kinetics energy : %f\n", individual_kin_energy);
        
        return velocities;
}

double compute_init_potential_energy(std::vector<Particle> particles, double L, double rc) {
    double potential_energy = 0.0;
    double rc6 = pow(rc, 6);
    double rc12 = rc6*rc6;
    int counter = 0;
    for (int i = 0; i < particles.size()-1; i++) {
        auto part = particles[i];
        double x = part.x, y = part.y;

        for (int j = i+1; j < particles.size(); j++) {
                auto part2 = particles[j];
                double d_x = x - part2.x;
                if (d_x > L/2) {
                    d_x -= (double) L;
                } else if (d_x <= -L/2) {
                    d_x +=  (double) L;
                }
                double d_y = y - part2.y;
                if (d_y > L/2) {
                    d_y -= (double)  L;
                } else if (d_y <= -L/2) {
                    d_y += (double)  L;
                }

                double r = sqrt(d_x*d_x + d_y*d_y);
                std::cout<< " r: " << r << endl;
                double r6 = pow(r, 6);
                double r12 = r6*r6;
                counter++;
                if (r <= rc) {
                    std::cout<< " Potential energy: " << potential_energy << endl;
                    //cout<< "value to add "<< 4 * (1.0/(r12) - 1.0/(r6)) - 4*((1.0/(rc12)) - (1.0/(rc6))) << endl;
                    potential_energy += 4 * (1.0/(r12) - 1.0/(r6)) - 4*((1.0/(rc12)) - (1.0/(rc6)));
                }
            
            
        }
    }
    
    std::cout << "unique interactions: " << counter << endl;
    return potential_energy;
}

double compute_init_kinetic_energy(std::vector<Particle> particles) {
    double kinetic_energy = 0.0;
    for (int i = 0; i < particles.size(); i++) {
        auto part = particles[i];
        kinetic_energy += particles[i].k_en;
    }
    return kinetic_energy;
}
int main(){
    

    int N;
    std::cout << "Insert N: \n";
    std::cin >> N;
    // Enable point size control in the vertex shader
    


    printf("Assignment 3\n");
    int L = 30;

    double m = 1.0;
    double rc = 2.5;
    double dt = 0.01;
    
    double thermostat_temp;
    std::cout << "Insert thermostat temperature: ";
    std::cin >> thermostat_temp;
    int N_X = L/rc;
    int N_Y = L/rc;
    printf("N_X = %d\n", N_X);
    thermostat_temp = 0.01;
    double x = 0;
    double y =0;
    double vx =0;
    double vy =0;
    double xcell =0;
    double ycell =0;
    std::vector<Particle> particle_list;



    printf("size vector: %lu", particle_list.size());
    double *velocities;
    velocities = return_velocities(0.01, N, m);
    
    for(int i=0; i<N; i++){
        std::cout << "velocity "<< i <<" = "<< velocities[i]<< std::endl;
    }
    std::vector<std::pair<double, double>> positions;
    for (int i=0;i<L;i++){
        for (int j=0;j<L;j++){
            positions.push_back({(double) i , (double) j});
        }
        
    };
    
    //for (int i=0;i<(positions.size());i++){
        //cout << "position "<< i << " = " << positions[i].first << ", "<<positions[i].second << endl;
        
    //};
    
    for (int i=0;i<N;i++){
        double alpha = 2.0 * 3.14 * rand();
        int random_indx = rand() % positions.size();
        double pos_x = positions[random_indx].first;
        double pos_y = positions[random_indx].second;
        positions.erase(positions.begin() + random_indx);
        positions.resize(positions.size());

        x = pos_x;
        y = pos_y;
        vx = velocities[i]*cos(alpha);
        vy = velocities[i]*sin(alpha);
        int xcell =   (int) x/  rc;
        int ycell = (int) y/  rc;
        cout << "xcell: " << xcell << " = " << (int) x<<" / "<< rc << endl;
        cout << "ycell: " << ycell << " = " << (int) y<<" / "<< rc << endl;
        particle_list.push_back(Particle(x, y, vx, vy,  m, L, xcell, ycell, thermostat_temp));
        
    };
    int counter =0;
    
    for (int i=0;i<particle_list.size();i++){
        cout << "particle "<< i << " position = " << particle_list[i].x << ", "<<particle_list[i].y << endl;
        //cout << "particle "<< i << " velocities = " << particle_list[i].vx << ", "<<particle_list[i].vy << endl;
        //cout << "particle "<< i << " cell = " << particle_list[i].i << ", "<<particle_list[i].j << endl;
    };
    double init_potential_energy = compute_init_potential_energy(particle_list, L, rc);
    cout << "Initial potential energy: " << init_potential_energy << endl;
    cout << "Initial kinetic energy: " << compute_init_kinetic_energy(particle_list) << endl;
    std::vector<std::vector<std::vector<Particle>>> cell_list(N_X, std::vector<std::vector<Particle>>(N_Y, std::vector<Particle>()));
    for (int i=0;i<particle_list.size();i++){
        cell_list[particle_list[i].i][particle_list[i].j].push_back(particle_list[i]);
    };

    
    vector<double>x_array = vector<double>(particle_list.size());
    vector<double>y_array = vector<double>(particle_list.size());
    for (int i=0;i<particle_list.size();i++){
        x_array[i] = particle_list[i].x;
        y_array[i] = particle_list[i].y;
    }
    

    
    zsock_t *socket = zsock_new_req("tcp://localhost:5555");
    zmsg_t *reply;
    while(1) {
        zmsg_t *message = zmsg_new();
        for (int i = 0; i < N_X; i++) {
            for (int j = 0; j < N_Y; j++) {
                for (auto it = cell_list[i][j].begin(); it != cell_list[i][j].end(); ) {
                    if (it->i == i && it->j == j) {
                        it = cell_list[i][j].erase(it);
                    }
                    else {
                        ++it;
                    }
                }
            }
        }
        
        int index =0;
        for (auto& part : particle_list) {
            part.update(dt, rc);
            int x_cell = static_cast<int>(part.x / rc);
            int y_cell = static_cast<int>(part.y / rc);
            cell_list[x_cell][y_cell].push_back(part);
            part.i = x_cell;
            part.j = y_cell;
            
            
        }
        
        for (int i=0;i<particle_list.size();i++){
            x_array[i] = particle_list[i].x;
            y_array[i] = particle_list[i].y;
        }
        zframe_t *frame1 = zframe_new(x_array.data(), sizeof(double) * x_array.size());
        zmsg_add(message, frame1);
        zframe_t *frame2 = zframe_new(y_array.data(), sizeof(double) * y_array.size());
        zmsg_add(message, frame2);
        zmsg_send(&message, socket);
        reply = zmsg_recv(socket);
        zframe_t *reply_frame = zmsg_pop(reply);
        std::cout << "Received reply: " << std::string((char *)(zframe_data(reply_frame)), zframe_size(reply_frame)) << std::endl;
        
        
        
    }
    zsock_destroy(&socket);
    zmsg_destroy(&reply);
    return 0;
}

