#include "./include/distrib2d.h"

double sqr(double x)
{
    return x*x;
}

std::vector<double> Distrib2D::read_electric_temp_coeffs(const std::string& path)
{
    std::ifstream in(path);
    std::vector<double> v{std::istream_iterator<double>{in}, std::istream_iterator<double>{}};
    
    return v;
}

double Distrib2D::eval_horners(double x) const
{
    //highest order coeff first
    double res = et_coeffs.front();

    for (int i = 1; i < et_coeffs.size(); i++)
    {
        res = res*x + et_coeffs[i];
    }

    return res;
}


// mV/m
double Distrib2D::electric_field(double t) const
{
    /*if (t <= 0.0)
    {
        return 0.0;
    }

    const double ramp = t*e_field_ramp_rate;
    if (ramp < e_field_peak)// rising edge
    {
        return ramp;
    }

    return std::max(0.0, e_field_peak - (ramp - e_field_peak));   // falling edge, min of 0
    */
    
    if (t <= p.th)
    {
        return 0.0;
    }
    if (t <= p.th + p.tau_h)
    {
        const double temp = p.Tc + rise*(t - p.th);
        const double e_val = eval_horners(std::log(temp));
        //const double e_val = eval_horners(temp);
        //std::print("{:.3f},{:.3f},{:.3f}\n", t, temp, e_val);
        return e_val;
    }
    if (t <= p.tc)
    {
        return p.Eh;
    }
    if (t <= p.tc + p.tau_c)
    {
        const double temp = p.Th + fall*(t - p.tc);
        const double e_val = eval_horners(std::log(temp));
        //const double e_val = eval_horners(temp);
        //std::print("{:.3f},{:.3f},{:.3f}\n", t, temp, e_val);
        return e_val;
    }
    if (t > p.tc + p.tau_c)
    {
        return 0.0;
    }

    std::print("{:.3f}\n", t);
    return -1.0;
}

//maybe precompute and interpolate
double Distrib2D::distribution_function(double v, double n, double t, int angle) const
{
    const double efield = electric_field(t);
    const double integral = ef.compute_integral(efield, angle);
    return n*ef.compute_dist(efield, angle, v)/integral;
}

/*
//check!
static double exobase_temp_par(const Distrib2D::SimTimeEvolutionParams& p, double t)
{
    if (t < p.th)
    {
        return p.Tc;
    }

    if (t < p.th + p.tau_h)
    {
        return p.Tc + (p.Th - p.Tc)*(t - p.th)/p.tau_h;
    }

    if (t < p.tc)
    {
        return p.Th;
    }

    if (t < p.tc + p.tau_c)
    {
        return p.Th + (p.Tc - p.Th)*(t - p.tc)/p.tau_c;
    }

    return p.Tc;
}

//check!
static double exobase_temp_perp(const Distrib2D::SimTimeEvolutionParams& p, double t)
{
    const double temp = exobase_temp_par(p, t);

    return (t >= p.th && t < p.tc + p.tau_c) ? Distrib2D::tpr*temp : temp;
}

//check!
double Distrib2D::distribution_function(double v, double n, double t, int angle) const
{
    const double temp = (angle == 0) ? exobase_temp_par(p, t) : exobase_temp_perp(p, t);
    const double b = mass/(2.0*kb*temp);

    return n*std::sqrt(b/pi)*std::exp(-b*v*v);
}*/

double Distrib2D::boundary_density(double t) const
{
    return 1e11;
    
    /*const double a = 1e12;

    if (t < 0.0)
    {
        return a;
    }

    return a + (t/1000.0)*a;*/

    //return 1e12;

    /*if (t < 0.0)
    {
        return 1e12;
    }
    return std::clamp(1e12*std::pow(0.5, t/200.0), 0.5e12, 1e12);*/
}

double Distrib2D::calc_dt(double vi, double vf, double t) const
{
    return (vf - vi)/acceleration;
}

double Distrib2D::compute_vi_par(double vf_par, double dz) const
{
    return std::sqrt(vf_par*vf_par - 2.0*acceleration*dz);
}

//optimize!
array2d<double> Distrib2D::get_f_vf_dist(double t, double dz) const
{
    array2d<double> result(nv, nv);

    #pragma omp parallel for
    for (int i = 0; i < nv; i++)
    {
        const double vf_par = vf_vec[i];
        const double vi_par = compute_vi_par(vf_par, dz);
        const double tau = t - calc_dt(vi_par, vf_par, t); //same particle
        const double bnd_density = boundary_density(tau);
        const double f_vf_par = distribution_function(vi_par, bnd_density, tau, 0);

        const double u = f_vf_par/bnd_density;

        for (int j = 0; j < nv; j++)
        {
            const double f_vf_perp = distribution_function(vf_vec[j], bnd_density, tau, 90);
            result[i, j] = u*f_vf_perp;
        }
    }

    return result;
}

//utility function
std::vector<double> Distrib2D::calc_times(double t1, double t2, double dt) const
{
    const int N = std::round((t2 - t1) / dt);
    std::vector<double> times(N);

    for (int i = 0; i < N; i++)
    {
        times[i] = t1 + i*dt;
    }

    return times;
}

/*Distrib2D::Moments Distrib2D::get_moments(double t1, double t2, double dt, double dz) const
{
    const int N = std::round((t2 - t1) / dt);
    Moments m(N);

    #pragma omp parallel for
    for (int n = 0; n < N; n++)
    {
        const double t = t1 + dt*n;

        double m0_par = 0.0;
        double m1_par = 0.0;
        double m2_par = 0.0;

        double m0_perp = 0.0;
        double m2_perp = 0.0;

        double integral0 = 0.0;
        double integral1 = 0.0;

        const double tau_perp = t - calc_dt(compute_vi_par(0.0, t, dz), 0.0, t);
        const double bnd_density_perp = boundary_density(tau_perp);

        for (int i = 0; i < nv; i++)
        {
            const double vf_par = vf_vec[i];
            const double vi_par = compute_vi_par(vf_par, dz);
            const double tau = t - calc_dt(vi_par, vf_par, t);
            const double bnd_density = boundary_density(tau);
            const double f_par = distribution_function(vi_par, bnd_density, tau, 0);
            const double f_perp = distribution_function(vf_vec[i], bnd_density_perp, tau_perp, 90); //vi = vf for perp
            const double wi = 1.0 - 0.5*((i == 0) + (i == nv - 1));

            const double v = vf_par;

            m0_par += wi*f_par;
            m1_par += wi*f_par*v;
            m2_par += wi*f_par*v*v;

            m0_perp += wi*f_perp;
            m2_perp += wi*f_perp*v*v;

            const double s = bnd_density / dv;

            integral0 += wi*s*f_par/bnd_density;
            integral1 += wi*s*f_par*vf_par/bnd_density;
        }

        const double density = dv*dv*integral0;
        const double avg_vel_par = integral1/integral0; //dv*dv*integral1 / density simplified

        m.density[n] = density;
        m.avg_par_velocity[n] = avg_vel_par;
        
        const double a = mass/kb;
        m.ion_temp_par[n] = a*(m2_par/m0_par - sqr(m1_par/m0_par));
        m.ion_temp_perp[n] = a*m2_perp/m0_perp;
    }

    return m;
}*/

/*
Distrib2D::Moments Distrib2D::get_moments(double t1, double t2, double dt, double dz) const
{
    const int N = std::round((t2 - t1) / dt);
    Moments m(N);

    std::vector<double> vi_par_vec(nv);

    for (int i = 0; i < nv; i++)
    {
        vi_par_vec[i] = compute_vi_par(vf_vec[i], dz);
    }

    #pragma omp parallel for schedule(dynamic)
    for (int n = 0; n < N; n++)
    {
        const double t = t1 + dt*n;

        double m0_par = 0.0;
        double m1_par = 0.0;
        double m2_par = 0.0;

        double m0_perp = 0.0;
        double m2_perp = 0.0;

        const double tau_perp = t - calc_dt(compute_vi_par(0.0, dz), 0.0, t);
        const double bnd_density_perp = boundary_density(tau_perp);

        for (int i = 0; i < nv; i++)
        {
            const double tau = t - calc_dt(vi_par_vec[i], vf_vec[i], t);
            const double bnd_density = boundary_density(tau);
            const double f_par = distribution_function(vi_par_vec[i], bnd_density, tau, 0);
            const double f_perp = distribution_function(vf_vec[i], bnd_density_perp, tau, 90); //vi = vf for perp

            const double wi = 1.0 - 0.5*((i == 0) + (i == nv - 1));

            const double v = vf_vec[i];

            m0_par += wi*f_par;
            m1_par += wi*f_par*v;
            m2_par += wi*f_par*v*v;

            m0_perp += wi*f_perp;
            m2_perp += wi*f_perp*v*v;
        }

        
        m.density[n] = dv*m0_par;
        m.avg_par_velocity[n] = m1_par/m0_par;

        const double a = mass/kb;
        m.ion_temp_par[n] = (m0_par > 0.0) ? a*(m2_par/m0_par - sqr(m1_par/m0_par)) : std::numeric_limits<double>::quiet_NaN();
        m.ion_temp_perp[n] = (m0_perp > 0.0) ? a*m2_perp/m0_perp : std::numeric_limits<double>::quiet_NaN();
    }

    return m;
}
*/


Distrib2D::Moments Distrib2D::get_moments(double t1, double t2, double dt, double dz) const
{
    constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
    const int N = std::round((t2 - t1) / dt);
    Moments m(N);

    std::vector<double> vi_par_vec(nv);

    for (int i = 0; i < nv; i++)
    {
        vi_par_vec[i] = compute_vi_par(vf_vec[i], dz);
    }

    #pragma omp parallel for schedule(dynamic)
    for (int n = 0; n < N; n++)
    {
        //std::print("{}/{}\n", n, N);
        const double t = t1 + dt*n;

        double m0 = 0.0;//sum f
        double m1_par = 0.0;//sum f*v_par
        double m2_par = 0.0;  //sum f*v_par^2
        double m1_perp = 0.0; //sum f*v_perp
        double m2_perp = 0.0; //sum f*v_perp^2

        double prev_efield = qnan;
        double s0 = 0.0;
        double s1 = 0.0;
        double s2 = 0.0;

        for (int i = 0; i < nv; i++)    
        {
            const double vf_par = vf_vec[i];
            const double tau = t - calc_dt(vi_par_vec[i], vf_par, t);
            const double bnd_density = boundary_density(tau);
            const double f_par = distribution_function(vi_par_vec[i], bnd_density, tau, 0);
            const double wi = 1.0 - 0.5*((i == 0) + (i == nv - 1));

            //perp part rides the same characteristic (tau set by parallel motion, v_perp
            //unchanged); unit-normalized so the product carries one factor of bnd_density
            const double efield = electric_field(tau);
            if (efield != prev_efield)
            {
                prev_efield = efield;
                const double integral = ef.compute_integral(efield, 90);

                s0 = 0.0;
                s1 = 0.0;
                s2 = 0.0;
                for (int j = 0; j < nv; j++)
                {
                    const double v_perp = vf_vec[j];
                    //const double f_perp = ef.compute_dist(efield, 90, v_perp)/integral;

                    const double f_perp = distribution_function(vf_vec[j], bnd_density, tau, 90);

                    const double wj = 1.0 - 0.5*((j == 0) + (j == nv - 1));

                    s0 += wj*f_perp;
                    s1 += wj*f_perp*v_perp;
                    s2 += wj*f_perp*v_perp*v_perp;
                }
            }

            const double g = wi*f_par;
            m0 += g*s0;
            m1_par += g*vf_par*s0;
            m2_par += g*vf_par*vf_par*s0;
            m1_perp += g*s1;
            m2_perp += g*s2;
        }

        //dv*dv cancels in every ratio, so only the density carries it
        const double avg_v_par = (m0 > 0.0) ? m1_par/m0 : qnan;
        const double avg_v_perp = (m0 > 0.0) ? m1_perp/m0 : qnan;

        m.density[n] = dv*dv*m0; //eq 14
        m.avg_par_velocity[n] = avg_v_par; //eq 15

        const double a = mass/kb;
        m.ion_temp_par[n] = (m0 > 0.0) ? a*(m2_par/m0 - sqr(avg_v_par)) : qnan;   //eq 16
        m.ion_temp_perp[n] = (m0 > 0.0) ? a*(m2_perp/m0 - sqr(avg_v_perp)) : qnan; //eq 17
    }

    return m;
}

void Distrib2D::test()
{
    for (double t = 0.0; t <= 4000.0; t += 0.5)
    {
        const double E = electric_field(t);
        const double inv = 1.0/ef.compute_integral(E, 90);
        double s0 = 0.0, s2 = 0.0;
        for (int j = 0; j < nv; j++)
        {
            const double v = vf_vec[j];
            const double f = ef.compute_dist(E, 90, v)*inv;
            const double wj = 1.0 - 0.5*((j == 0) + (j == nv - 1));
            s0 += wj*f;
            s2 += wj*f*v*v;
        }
        std::print("{:.0f},{:.3f},{:.1f}\n", t, E, (mass/kb)*s2/s0);
    }
}

void Distrib2D::boundary_sweep()
{
    const double a = mass/kb;
    for (double t = 0.0; t <= 4000.0; t += 4.0)
    {
        const double E = electric_field(t);
        double s0[2] = {0.0, 0.0}, s2[2] = {0.0, 0.0};

        for (int k = 0; k < 2; k++)
        {
            const int ang = 90*k;
            const double inv = 1.0/ef.compute_integral(E, ang);
            for (int j = 0; j < nv; j++)
            {
                const double v = vf_vec[j];
                const double f = ef.compute_dist(E, ang, v)*inv;
                const double wj = 1.0 - 0.5*((j == 0) + (j == nv - 1));
                s0[k] += wj*f;
                s2[k] += wj*f*v*v;
            }
        }

        std::print("{:.0f},{:.4f},{:.1f},{:.1f},{:.4f}\n",
                   t, E, a*s2[0]/s0[0], a*s2[1]/s0[1], dv*s0[1]);
    }
}