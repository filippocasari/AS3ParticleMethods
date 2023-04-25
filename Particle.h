#include <cmath>
class Particle
{
public:
    double x;
    double y;
    double thermostat_temp;
    double vx;
    double vy;
    double m;
    int L;
    double a_x = 0.;
    double a_y = 0.;
    double F_x = 0.;
    double F_y = 0.;
    int i;
    int j;
    double k_en;
    Particle(double x, double y, double vx, double vy, double m, int L, int i, int j, double thermostat_temp)
    {
        this->x = x;
        this->y = y;
        this->vx = vx;
        this->vy = vy;
        this->m = m;
        this->L = L;
        this->i = i;
        this->j = j;
        this->thermostat_temp = thermostat_temp;
        this->compute_k_energy();
    }
    void update(double dt, double rc)
    {
        this->y += (this->vy * dt) + this->a_y * dt * dt * 0.5;
        this->x += (this->vx * dt) + this->a_x * dt * dt * 0.5;
        this->x = fmod(this->x, this->L);
        this->y = fmod(this->y, this->L);
        this->i = int(this->x / rc);
        this->j = int(this->y / rc);
    }
    void update_vel_acc(double dt, double F_x, double F_y)
    {

        this->F_x = F_x;
        this->F_y = F_y;
        double new_acc_x = this->F_x / this->m;
        double new_acc_y = this->F_y / this->m;
        this->vx += (this->a_x + new_acc_x) * dt * 0.5;
        this->vy += (this->a_y + new_acc_y) * dt * 0.5;
        this->a_x = new_acc_x;
        this->a_y = new_acc_y;
        this->compute_k_energy();
    }
    void compute_k_energy()
    {
        this->k_en = 0.5 * this->m * (pow(this->vx, 2) + pow(this->vy, 2));
    }
    bool operator==(const Particle& other) const {
        return x == other.x && y == other.y && vx == other.vx && vy == other.vy && m == other.m && L == other.L && i == other.i && j == other.j;
    }
    bool operator!=(const Particle& other) const {
        return !(*this == other);
    }
};