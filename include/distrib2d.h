#pragma once
#include <cmath>
#include <vector>
#include <print>
#include <numbers>
#include <climits>

#include "mdarray.h"
#include "ElectricField.h"

class Distrib2D
{
public:
    //moments at time steps
    struct Moments
    {
        Moments(int N) : density(N), avg_par_velocity(N),
                         ion_temp_par(N) , ion_temp_perp(N){}

        std::vector<double> density;
        std::vector<double> avg_par_velocity;
        std::vector<double> ion_temp_par;
        std::vector<double> ion_temp_perp;
    };

    //TODO: add Ec
    struct SimTimeEvolutionParams
    {
        double Tc, Th, tc, th, tau_c, tau_h, Eh;
        std::string et_coeff_path;
    };

    Distrib2D(double v_min, double v_max, double v_step, double m, std::string coeff_path,
        const Distrib2D::SimTimeEvolutionParams& params) : 
        vf_min(v_min), vf_max(v_max), dv(v_step), mass(m), ef(coeff_path), p(params)
    {
        nv = std::round((vf_max - vf_min) / dv) + 1;

        //populate vf arrays
        vf_vec.resize(nv);
        for (int i = 0; i < nv; i++)
        {
            vf_vec[i] = vf_min + dv*i;
        }

        et_coeffs = read_electric_temp_coeffs(p.et_coeff_path);
        rise = (p.Th - p.Tc)/p.tau_h;
        fall = (p.Tc - p.Th)/p.tau_c;
        p.Eh = eval_horners(std::log(p.Th));
        //p.Eh = eval_horners(p.Th);

        std::cout << p.Eh << '\n';

        test();
        boundary_sweep();
    }

    array2d<double> get_f_vf_dist(double t, double dz) const;
    Moments get_moments(double t1, double t2, double dt, double dz) const;
    std::vector<double> calc_times(double t1, double t2, double dt) const;
    double electric_field(double t) const; // [mV/m] only here for testing
    void test();
    void boundary_sweep();

    //change
    static constexpr double tpr = 1.2;

private:
    double vf_min, vf_max, dv, mass;
    std::vector<double> vf_vec;
    int nv;
    ElectricField ef;
    std::vector<double> et_coeffs;
    SimTimeEvolutionParams p;
    double rise, fall;

    static constexpr double pi = std::numbers::pi_v<double>;
    static constexpr double kb = 1.380649e-23; //https://physics.nist.gov/cgi-bin/cuu/Value?k
    static constexpr double acceleration = -4.905; //m s^-2
    static constexpr double e_field_ramp_rate = 1.0; //1 mV/m per second
    static constexpr double e_field_peak = 100.0; //  [mV/m]

    double eval_horners(double x) const;
    std::vector<double> read_electric_temp_coeffs(const std::string& path);
    
    double boundary_density(double t) const; 

    double distribution_function(double v, double n, double t, int angle) const;

    double calc_dt(double vi, double vf, double t) const;
    double compute_vi_par(double vf_par, double dz) const;
};