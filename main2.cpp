#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <print>
#include <filesystem>

#include <fenv.h>

#include "./include/mdarray.h"
#include "./include/distrib2d.h"
#include "./include/npy.hpp" // https://github.com/llohse/libnpy


void write_vec1d_npy(const std::string& filename, const std::vector<double>& vec)
{
    npy::npy_data_ptr<double> d;
    d.data_ptr = vec.data();
    d.shape = {vec.size()};
    d.fortran_order = false;

    npy::write_npy(filename, d);
}

void write_array2d_npy(const std::string& filename, const array2d<double>& arr)
{
    npy::npy_data_ptr<double> d;
    d.data_ptr = arr.data();
    d.shape = {arr.dim0(), arr.dim1()};
    d.fortran_order = false;

    npy::write_npy(filename, d);
}

Distrib2D::SimTimeEvolutionParams process_args(char* argv[])
{
    Distrib2D::SimTimeEvolutionParams params;
    params.et_coeff_path = argv[2];
    params.Tc = std::stod(argv[3]);
    params.Th = std::stod(argv[4]);
    params.tc = std::stod(argv[5]);
    params.th = std::stod(argv[6]);
    params.tau_c = std::stod(argv[7]);
    params.tau_h = std::stod(argv[8]);

    return params;
}

//args: coeff dir, E(T) coeff path, Tc, Th, tc, th, tau_c, tau_h, z1, z2, dz
int main(int argc, char* argv[])
{
    feenableexcept(FE_INVALID);

    if (argc < 2)
    {
        std::cerr << "invalid arg count\n";
        return 1; 
    }

    const std::string dir(argv[1]); //brace init? ill look into it

    if (dir == "-h" || dir == "--help")
    {
        std::cout << "args: coeff dir, E(T) coeff path, Tc, Th, tc, th, tau_c, tau_h, z1, z2, dz[km]\n";
        return 0;
    }

    if (argc != 12)
    {
        std::cerr << "invalid arg count\n";
        return 1; 
    }

    Distrib2D::SimTimeEvolutionParams params = process_args(argv);

    const std::string moment_dir = "./moment_output/";
    const std::string grid_dir = "./output_2d/";

    std::filesystem::path md{moment_dir};
    std::filesystem::path gd{grid_dir};

    const bool created_md = std::filesystem::create_directories(md);
    const bool created_gd = std::filesystem::create_directories(gd);

    if (created_md)
    {
        std::print("created moment output directory\n");
    }

    if (created_gd)
    {
        std::print("created grid output directory\n");
    }

    const double z1 = 1000.0*std::stod(argv[9]);
    const double z2 = 1000.0*std::stod(argv[10]);
    const double dz = 1000.0*std::stod(argv[11]);
    const int nz = std::round((z2 - z1)/dz) + 1;
    const double mass = 2.66e-26;
    
    const double t1 = 0.0;
    const double t2 = 4000.0;
    const double dt = 4.0;

    const double vmin = -10500.0;
    const double vmax = -vmin;
    const double dv = 10.0;

    Distrib2D dist{vmin, vmax, dv, mass, dir, params};
    std::vector<double> times = dist.calc_times(t1, t2, dt);

    double p1 = 1000.0;
    double p2 = 1000.0;

    for (int i = 0; i < nz; i++)
    {
        const auto start = std::chrono::steady_clock::now();

        const double z = z1 + i*dz;

        Distrib2D::Moments m = dist.get_moments(t1, t2, dt, z);

        const std::string dz_str = std::to_string(static_cast<int>(z/1000)) + "km";
        write_vec1d_npy(moment_dir + dz_str + "_time.npy", times);
        write_vec1d_npy(moment_dir + dz_str + "_density.npy", m.density);
        write_vec1d_npy(moment_dir + dz_str + "_avg_par_velocity.npy", m.avg_par_velocity);
        write_vec1d_npy(moment_dir + dz_str + "_ion_temp_par.npy", m.ion_temp_par);
        write_vec1d_npy(moment_dir + dz_str + "_ion_temp_perp.npy", m.ion_temp_perp);

        const auto end = std::chrono::steady_clock::now();
        const int ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        const double percent = 100.0*(i + 1.0)/static_cast<double>(nz);
        const double est_rem = (nz - i)*(ms + p1 + p2)/(3.0*1000.0);
        p1 = ms;
        p2 = p1;


        std::print("dz = {:.1f}[km]\t{:.2f}[s] \t estimated time remaining: {:.1f}[mins] \t {}/{} iters \t {:.1f}%\n", z/1000.0, ms/1000.0, est_rem/60.0, i + 1, nz, percent);
    }



    return 0; 


}



    //
    //int ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    //std::print("dz = {:.1f}[km], {}[ms]        \n", dz/1000.0, ms);
    /*std::vector<double> e_vals(times.size());
    for (int i = 0; i < e_vals.size(); i++)
    {
        e_vals[i] = dist.electric_field(times[i]);
    }

    write_vec1d_npy("./electric_field_vals.npy", e_vals);*/