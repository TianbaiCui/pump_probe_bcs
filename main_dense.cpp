//
//  main.cpp
//  THz-Pump-Probe
//
//  Created by Tianbai Cui on 7/11/18.
//  Copyright © 2018 Tianbai Cui. All rights reserved.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <fstream>
#include <string>
#include <boost/math/special_functions/ellint_1.hpp>

#include "nr3.h"
#include "Fortran.h"
#include "Matrix.h"
#include "interp_1d.h"
#include "quadrature.h"
#include "romberg.h"
#include "stepper.h"
#include "odeint.h"
#include "stepperdopr5.h"
#include "stepperdopr853.h"
#include "roots.h"
#include "params.h"

using namespace std;
const double pi = 3.14159265358979;

/********************************************************************
 Define functions for calculating the thermalized gap and temperature
 ********************************************************************/

double fu_x (double x, double u) {
    return 2. * cosh(2. * x) * (1. - tanh(u * cosh(x)));
}

double log_gap (double x, double u) {
    return 2./(sqrt(x*x + 4.*u*u) * (1. + exp(sqrt(x*x + 4.*u*u))));
}

double Fu (double u) {
    auto fu_x2 = [&u](double x) {
        return fu_x(x, u);
    };
    double Int_x = qromb(fu_x2, 0., 4 * omega_D / Delta0);
    return Int_x;
}

double gap (double y) {
    auto log_gap_u = [&y](double x) {
        return log_gap(x, y);
    };
    double Int_x_gap = qromb(log_gap_u, 0. , 4 * omega_D / Delta0);
    return exp(-Int_x_gap);
}

double solve_u_E (double u, double En) {
    return pow(gap(u),2) * (1. - Fu(u)) - 1. + 2. * En;
}

double solve_0_E (double T, double En) {
    auto integrand = [&T] (double x) {
        return sqrt(x*x+1.*1.) - x * tanh(x / (2. * T));
    };
    double Int_x_0 = qromb(integrand, 0., omega_D / Delta0);
    return 2. * Int_x_0 - 1/(V0*rho0) - En;  //log(2.*omega_D / Delta0)
}

double En_u (double u) {
    return (- pow(gap(u),2) * (1. - Fu(u)) + 1.)/2.;
}

double calculate_T (double en) {
    double sol_u_temp, Delta_eq_temp, T_eq_temp;
    double en_Tianbai = en;

    if (en_Tianbai < 1.0e-8) {
        sol_u_temp = 1/0.04;
        Delta_eq_temp = 1. * Delta0;
        T_eq_temp = 0.04 * Delta0;
    }
    else if (en_Tianbai > 1.5574) {
        sol_u_temp = 0;
        Delta_eq_temp = 0.;
        auto solve_0 = [&en_Tianbai](double T) {
            return solve_0_E(T, en_Tianbai);
        };
        T_eq_temp = zbrent(solve_0, 0.566933, 10, 1.0e-5) * Delta0;
    }
    else {
        auto solve_u = [&en_Tianbai](double u) {
            return solve_u_E(u,en_Tianbai);
        };
        sol_u_temp = zbrent(solve_u, 0.001, 10., 1.0e-5);
        Delta_eq_temp = gap(sol_u_temp) * Delta0;
        T_eq_temp = Delta_eq_temp / (2. * sol_u_temp);
    }
    return T_eq_temp;
}

/********************************************************************************
 Define functions to calculate the equilibrium gap value using energy integration
 ********************************************************************************/
double gap_eq (double delta) {
    double output = 0.;
    double DOS;
    for (int idx = 0; idx < N; idx++) {
        double xi = - omega_D + 2*omega_D/(static_cast<double>(N-1))*static_cast<double>(idx);
        DOS = 1/(2*pi*pi*J)*boost::math::ellint_1(sqrt(1-pow((xi + mu)/(2*J), 2)/4));
        output += (V0*DOS*2.*omega_D/N)/(2*sqrt(xi*xi + delta*delta))*tanh(sqrt(pow(delta,2) + pow(xi,2))/(2 * 0.04 * delta));
    }
    return output-1.;
}

/********************************************************************************
 Calculate the equilibrium gap value using momentum summation
 ********************************************************************************/
double gap_eq_sum_k (double delta){
    double output = 0.;
    for (int idx = 0; idx < N; idx++) {
        double kx = pi/(static_cast<double>(N-1))*static_cast<double>(idx);
        for (int idy = 0; idy < N; idy++) {
            double ky = pi/(static_cast<double>(N-1))*static_cast<double>(idy);
            double xi = - 2*J * (cos(kx)+cos(ky)) - mu;
            if (xi >= -omega_D && xi <= omega_D){
                output += 1/(2*sqrt(xi*xi+delta*delta))*(V0/(N*N))*tanh(sqrt(pow(delta,2) + pow(xi,2))/(2 * 0.04 * delta));
            }
        }
    }
    return output - 1.;
}

/********************************************************************************
 Calculate the equilibrium gap value using Peter's momentum discretization
 ********************************************************************************/
double gap_eq_en_k (double delta){
    double output = 0.;
    Int N_k = L_energy * L_kx;
    VecDoub kstates(3*L_energy*L_kx);
    for(int idx_en = 0; idx_en < L_energy; idx_en++){
        //double en = mu-omega_D+2*omega_D*idx_en/static_cast<double>(L_energy-1);
        double en = mu-omega_D+2*J*(cos(A0)-1)+((1/cos(A0)+1)*omega_D+2*J*(1/cos(A0)-cos(A0))+mu*(1/cos(A0)-1))*idx_en/static_cast<double>(L_energy-1);
        double kx_min, kx_max, kx, ky, ky_max;
        if (en < 0.){
            kx_min = 0.;
            kx_max = acos(-1 + abs(en)/(2.*J));
        }
        else {
            kx_min = acos(1-en/(2.*J));
            kx_max = pi;
        }
        for (int idx_kx = 0.; idx_kx < L_kx/2; idx_kx++){
            kx = kx_min + (7.*kx_max/10. - kx_min)*idx_kx/static_cast<double>(L_kx/2-1);
            if (en < 0.){
                ky = acos(abs(en)/(2.*J) - cos(kx));
            }
            else{
                ky = acos(-en/(2.*J) - cos(kx));
            }
            kstates[3*idx_en*L_kx + 3*idx_kx] = en;
            kstates[3*idx_en*L_kx + 3*idx_kx + 1] = kx;
            kstates[3*idx_en*L_kx + 3*idx_kx + 2] = ky;
        }
        ky_max = ky;
        for (int idx_kx = L_kx/2; idx_kx < L_kx; idx_kx++){
            ky = kx_min + (ky_max - kx_min)*(idx_kx - L_kx/2)/static_cast<double>(L_kx/2);
            if (en < 0.) {
                kx = acos(abs(en)/(2.*J) - cos(ky));
            }
            else {
                kx = acos(-en/(2.*J) - cos(ky));
            }
            kstates[3*idx_en*L_kx + 3*idx_kx] = en;
            kstates[3*idx_en*L_kx + 3*idx_kx + 1] = kx;
            kstates[3*idx_en*L_kx + 3*idx_kx + 2] = ky;
        }
    }// end momentum grid
    for (int idx = 0; idx < N; idx++) {
        double xi = kstates[3*idx];
        if (xi>=-omega_D && xi<=omega_D){
            output += 1/(2*sqrt(xi*xi+delta*delta))*(V0/N)*tanh(sqrt(pow(delta,2) + pow(xi,2))/(2 * 0.04 * delta));//*2*omega_D/L_energy;
            //output += V0/(2*sqrt(xi*xi + delta*delta))*tanh(sqrt(pow(delta,2) + pow(xi,2))/(2 * 0.04 * delta))/N_k;
        }
    }
    return output - 1.;
}

/********************************************************************************
 Calculate the equilibrium gap value using a trianglar momentum discretization
 ********************************************************************************/
double gap_eq_sum_kxy (double delta, int L_kxy) {
    double output = 0.;
    //double A0 = 0.;
    for (int idx = 0; idx < L_kxy; idx++) {
        double kxy = acos((omega_D-mu)/(2*J)-cos(A0)) + (2*acos((omega_D+mu)/(-4*J*cos(A0))-0.5*(1/cos(A0)-1))-acos((omega_D-mu)/(2*J)-cos(A0)))/(static_cast<double>(L_kxy-1))*static_cast<double>(idx);
        //double kxy = acos((omega_D-mu)/(2*J)) + (sqrt(2)*acos(-(omega_D+mu)/(4*J))-acos((omega_D-mu)/(2*J)))/(static_cast<double>(N-1))*static_cast<double>(idx);
        for (int idy = 0; idy < L_kx; idy++) {
            double ky = 2*acos((omega_D+mu)/(-4*J*cos(A0))-0.5*(1/cos(A0)-1))/(static_cast<double>(L_kx-1))*static_cast<double>(idy);
            double xi = - 2*J * (cos(kxy-ky) + cos(ky)) - mu;
            if (xi >= -omega_D && xi <= omega_D && kxy-ky >= 0){
                output += 1./(2.*sqrt(xi*xi+delta*delta))*(V0/(L_kxy*L_kx))*tanh(sqrt(pow(delta,2) + pow(xi,2))/(2 * 0.04 * delta));
            }
        }
    }
    return output - 1.;
}

/*******************************************************************
 Define the structures for ODE solver
 ********************************************************************/
struct rhs_bcs_k_space {
    int N;
    double V0;
    double J;
    double mu;
    double omega_D;
    double Delta0;
    double A0x;
    double A0y;
    double tau;
    double sigma;
    double omega_pump;
    double T1, T2;
    VecDoub kstates;
    double ps_ratio;

    rhs_bcs_k_space(int NN, double VV0, double JJ, double mmu, double oomega_D, double DDelta0, double AA0x, double AA0y, double ttau, double ssigma, double oomega_pump, double TT1, double TT2, VecDoub kkstates, double pps_ratio) : N(NN), V0(VV0), J(JJ), mu(mmu), omega_D(oomega_D), Delta0(DDelta0), A0x(AA0x), A0y(AA0y), tau(ttau), sigma(ssigma), omega_pump(oomega_pump), T1(TT1), T2(TT2), kstates(kkstates), ps_ratio(pps_ratio){}
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        // define vector potential of laser field at time x
        double Ax=0.;
        double Ay=0.;
        double eta = 1.;
        //double A = 0.;
        if ((x > 0.) && (x < 2.*tau)) {
            Ax = A0x*exp(-pow(x-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*x);
            Ay = A0y*exp(-pow(x-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*x);
            //A = A0*exp(-pow(x-tau/2., 2)/(2.*sigma*sigma) ) * cos(omega_pump*x);
        }
        // first calculate the gap and energy for the current spin state y = S^\alpha_k
        double DeltaRe = 0.;
        double DeltaIm = 0.;
        double energy = 0.;
        for (int idx = 0; idx < N; idx++) {
            double kx = kstates[3*idx + 1];
            double ky = kstates[3*idx + 2];
            //double xi = - 2*J * (cos(kx) * cos(Ax) + cos(ky) * cos(Ay)) - mu;
            //if (xi >= -omega_D && xi <= omega_D && kx >=0){

                //energy += (ps_ratio/N) * (2. * xi * y[3*idx + 2] + xi);
            //}
            if (kstates[3*idx] >= -omega_D && kstates[3*idx] <= omega_D && kx>=0){
                DeltaRe += ps_ratio*y[3*idx]*(V0/N);
                DeltaIm += ps_ratio*y[3*idx + 1]*(-V0/N);
                energy += (ps_ratio/N) * (2. * kstates[3*idx] * y[3*idx + 2] + sqrt(kstates[3*idx]*kstates[3*idx] + Delta0*Delta0));
            }
            //if (kstates[3*idx] >= -omega_D && kstates[3*idx] <= omega_D && kx>=0){
            //    energy += (ps_ratio/N) * (2. * xi * y[3*idx + 2] + xi - kstates[3*idx] + sqrt(kstates[3*idx]*kstates[3*idx] + Delta0*Delta0));
            //}
        }
        energy += -1./V0*(DeltaRe*DeltaRe + DeltaIm*DeltaIm + Delta0*Delta0);
        energy = energy/(Delta0*Delta0);
        // determine T_eq corresponding to this energy
        double sol_u, Delta_eq, T_eq;
        double energy_TB = energy/rho0;
        if (energy_TB < 1.0e-8) {
            sol_u = 1/0.04;
            Delta_eq = 1. * Delta0;
            T_eq = 0.04 * Delta0;
        }
        else if (energy_TB > 1.5574) {
            sol_u = 0;
            Delta_eq = 0.;
            auto solve_0 = [&energy_TB](double T) {
                return solve_0_E(T, energy_TB);
            };
            T_eq = zbrent(solve_0, 0.56, 20, 1.0e-5) * Delta0;//0.566933
        }
        else {
            auto solve_u = [&energy_TB](double u) {
                return solve_u_E(u,energy_TB);
            };
            sol_u = zbrent(solve_u, 0.001, 10., 1.0e-5);
            Delta_eq = gap(sol_u) * Delta0;
            T_eq = Delta_eq / (2. * sol_u);
        }
        //cout << "energy/Delta0^2 = " << energy << ", T_eq = " << T_eq <<  ", T/Tc = " << T_eq/(Delta0*0.566933) << ", Delta_eq = " << Delta_eq << ", DeltaRe = " << DeltaRe << ", DeltaRe/Delta0 = " << DeltaRe/Delta0 << endl;

        // define the rhs of the ODE
        for (int idx = 0; idx < N; idx++) {
            double kx = kstates[3*idx+1];
            double ky = kstates[3*idx+2];
            double en = kstates[3*idx];
            //double xi = -2*J * (cos(kx) * cos(Ax) + cos(ky)*cos(Ay)) - mu;
            double xi = -mu;
            double D_Re = 0.;
            double D_Im = 0.;
            if (kx>=0){
                xi = -2*J*(cos(kx)*cos(Ax)+cos(ky)*cos(Ay))-mu;
                D_Re = DeltaRe;
                D_Im = DeltaIm;
            }
            double E_bog = sqrt(en*en + Delta_eq*Delta_eq);
            double sin_eps, cos_eps, sin_phi, cos_phi;
            if (Delta_eq > 1.e-6 && kx>=0) {
                sin_eps = Delta_eq/E_bog;
                cos_eps = -en/E_bog;
                cos_phi = DeltaRe / sqrt(DeltaRe*DeltaRe+DeltaIm*DeltaIm);
                sin_phi = DeltaIm / sqrt(DeltaRe*DeltaRe+DeltaIm*DeltaIm);
            }
            else if (en < 0.) {
                sin_eps = 0.;
                cos_eps = 1.;
                cos_phi = 1.;
                sin_phi = 0.;
            }
            else {
                sin_eps = 0.;
                cos_eps = -1.;
                cos_phi = 1.;
                sin_phi = 0.;
            }
            if ((T1 > 1e-5) && (T2 > 1e-5)) {
                dydx[3*idx] = 2.*(-xi * y[3*idx + 1] + D_Im * y[3*idx + 2]) - sin_eps*cos_phi/T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq))) - sin_phi/T2 * (sin_phi * y[3*idx] + cos_phi * y[3*idx+1])- cos_eps*cos_phi/T2 * (cos_eps*cos_phi * y[3*idx] - cos_eps*sin_phi * y[3*idx+1] - sin_eps * y[3*idx + 2]);

                dydx[3*idx + 1] = 2.*(D_Re*y[3*idx + 2] + xi * y[3*idx]) + sin_phi*sin_eps/T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq)))- cos_phi/T2 * (sin_phi * y[3*idx] + cos_phi*y[3*idx + 1]) - cos_eps*sin_phi/T2 * (-cos_eps*cos_phi*y[3*idx] + cos_eps*sin_phi*y[3*idx+1] + sin_eps*y[3*idx+2]);

                dydx[3*idx + 2] = -2.*(D_Im * y[3*idx] + D_Re * y[3*idx + 1]) - cos_eps/T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq))) - sin_eps/T2 * (-cos_eps*cos_phi * y[3*idx] + cos_eps*sin_phi * y[3*idx+1] + sin_eps * y[3*idx + 2]);
            }
            else if ((T1 < 1e-5) && (T2 > 1e-5)) {
              dydx[3*idx] = 2.*(-xi * y[3*idx + 1] + D_Im * y[3*idx + 2]) - sin_phi/T2 * (sin_phi * y[3*idx] + cos_phi * y[3*idx+1])- cos_eps*cos_phi/T2 * (cos_eps*cos_phi * y[3*idx] - cos_eps*sin_phi * y[3*idx+1] - sin_eps * y[3*idx + 2]);

              dydx[3*idx + 1] = 2.*(D_Re*y[3*idx + 2] + xi * y[3*idx]) - cos_phi/T2 * (sin_phi * y[3*idx] + cos_phi*y[3*idx + 1]) - cos_eps*sin_phi/T2 * (-cos_eps*cos_phi*y[3*idx] + cos_eps*sin_phi*y[3*idx+1] + sin_eps*y[3*idx+2]);

              dydx[3*idx + 2] = -2.*(D_Im * y[3*idx] + D_Re * y[3*idx + 1]) - sin_eps/T2 * (-cos_eps*cos_phi * y[3*idx] + cos_eps*sin_phi * y[3*idx+1] + sin_eps * y[3*idx + 2]);
            }
            else if ((T1 > 1e-5) && (T2 < 1e-5)) {
              dydx[3*idx] = 2.*(-xi * y[3*idx + 1] + D_Im * y[3*idx + 2]) - sin_eps*cos_phi/T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq)));

              dydx[3*idx + 1] = 2.*(D_Re*y[3*idx + 2] + xi * y[3*idx]) + sin_phi*sin_eps/T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq)));

              dydx[3*idx + 2] = -2.*(D_Im * y[3*idx] + D_Re * y[3*idx + 1]) - cos_eps/T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq)));
            }
            else {
                dydx[3*idx] = 2.*(-xi * y[3*idx +1] + D_Im * y[3*idx + 2]) ;

                dydx[3*idx + 1] = 2.*(D_Re*y[3*idx +2] + xi * y[3*idx]) ;

                dydx[3*idx + 2] = 2.*(-D_Im * y[3*idx] - D_Re * y[3*idx + 1]);
            }
        } // closes the for-loop
    } // closes void operator ()
};

struct rhs_bcs_k_space2 {
    int N;
    double V0;
    double J;
    double mu;
    double omega_D;
    double Delta0;
    double T_eq_final, Delta_eq_final;
    double DeltaRe_final, DeltaIm_final;
    double T1, T2;
    VecDoub kstates;
    double ps_ratio;

    rhs_bcs_k_space2(int NN, double VV0, double JJ, double mmu, double oomega_D, double DDelta0, double TT_eq_final, double DDelta_eq_final, double DDeltaRe_final, double DDeltaIm_final, double TT1, double TT2, VecDoub kkstates, double pps_ratio) : N(NN), V0(VV0), J(JJ), mu(mmu), omega_D(oomega_D), Delta0(DDelta0), T_eq_final(TT_eq_final), Delta_eq_final(DDelta_eq_final), DeltaRe_final(DDeltaRe_final), DeltaIm_final(DDeltaIm_final), T1(TT1), T2(TT2), kstates(kkstates), ps_ratio(pps_ratio){}
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx) {
        // first calculate the gap and energy for the current spin state y = S^\alpha_k
        double DeltaRe = 0.;
        double DeltaIm = 0.;
        //double energy = 0.;
        for (int idx = 0; idx < N; idx++) {
            double xi = kstates[3*idx];
            if (xi>=-omega_D && xi<=omega_D && kstates[3*idx+1]>=0){
                DeltaRe += ps_ratio*y[3*idx]*(V0/N);
                DeltaIm += ps_ratio*y[3*idx+1]*(-V0/N);
                //energy += (ps_ratio/N) * (2.* xi * y[3*idx + 2] + xi - kstates[3*idx] + sqrt(kstates[3*idx]*kstates[3*idx] + Delta0*Delta0));
            }
            //energy += (1/(N*N)) * (2.*xi*y[3*idx*N + 3*idy + 2] + sqrt(xi*xi + Delta0*Delta0));
            //energy += (1./N) * (2.*xi*y[3*idx+2] + 2.* DeltaRe * y[3*idx] - 2.* DeltaIm * y[3*idx+1] + sqrt(xi*xi + Delta0*Delta0));
        }
        // define the rhs of the ODE
        for (int idx = 0; idx < N; idx++) {
            double xi = kstates[3*idx];
            double D_Re = DeltaRe;
            double D_Im = DeltaIm;
            if (kstates[3*idx+1]<0){
                xi = -mu;
                D_Re = 0.;
                D_Im = 0.;
            }
            double E_bog = sqrt(xi*xi + Delta_eq_final*Delta_eq_final);
            double sin_eps, cos_eps, sin_phi, cos_phi;
            if (Delta_eq_final > 1.e-6 && kstates[3*idx+1]>=0) {
                sin_eps = Delta_eq_final/E_bog;
                cos_eps = -xi/E_bog;
                cos_phi = DeltaRe/sqrt(DeltaRe*DeltaRe + DeltaIm*DeltaIm);
                sin_phi = DeltaIm/sqrt(DeltaRe*DeltaRe + DeltaIm*DeltaIm);
            }
            else if (xi < 0.) {
                sin_eps = 0.;
                cos_eps = 1.;
                cos_phi = 1.;
                sin_phi = 0.;
            }
            else {
                sin_eps = 0.;
                cos_eps = -1.;
                cos_phi = 1.;
                sin_phi = 0.;
            }
            if ((T1 > 1e-5) && (T2 > 1e-5)) {
                dydx[3*idx] = 2.*(-xi * y[3*idx + 1] + D_Im * y[3*idx + 2]) - sin_eps * cos_phi / T1 * (sin_eps * cos_phi * y[3*idx] - sin_eps * sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq_final))) - sin_phi / T2 * (sin_phi * y[3*idx] + cos_phi * y[3*idx+1]) - cos_eps * cos_phi / T2 * (cos_eps * cos_phi * y[3*idx] - cos_eps * sin_phi * y[3*idx+1] - sin_eps * y[3*idx + 2]);

                dydx[3*idx + 1] = 2.*(D_Re*y[3*idx + 2] + xi * y[3*idx]) + sin_phi * sin_eps / T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq_final)))- cos_phi / T2 * (sin_phi * y[3*idx] + cos_phi * y[3*idx + 1]) - cos_eps * sin_phi/T2 * (-cos_eps * cos_phi * y[3*idx] + cos_eps * sin_phi * y[3*idx+1] + sin_eps * y[3*idx+2]);

                dydx[3*idx + 2] = -2.*(D_Im * y[3*idx] + D_Re * y[3*idx + 1]) - cos_eps / T1 * (sin_eps * cos_phi * y[3*idx]- sin_eps * sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq_final))) - sin_eps / T2 * (-cos_eps * cos_phi * y[3*idx] + cos_eps * sin_phi * y[3*idx+1] + sin_eps * y[3*idx + 2]);
            }
            else if ((T1 < 1e-5) && (T2 > 1e-5)) {
              dydx[3*idx] = 2.*(-xi * y[3*idx + 1] + D_Im * y[3*idx + 2]) - sin_phi / T2 * (sin_phi * y[3*idx] + cos_phi * y[3*idx+1]) - cos_eps * cos_phi / T2 * (cos_eps * cos_phi * y[3*idx] - cos_eps * sin_phi * y[3*idx+1] - sin_eps * y[3*idx + 2]);

              dydx[3*idx + 1] = 2.*(D_Re*y[3*idx + 2] + xi * y[3*idx]) - cos_phi / T2 * (sin_phi * y[3*idx] + cos_phi * y[3*idx + 1]) - cos_eps * sin_phi/T2 * (-cos_eps * cos_phi * y[3*idx] + cos_eps * sin_phi * y[3*idx+1] + sin_eps * y[3*idx+2]);

              dydx[3*idx + 2] = -2.*(D_Im * y[3*idx] + D_Re * y[3*idx + 1]) - sin_eps / T2 * (-cos_eps * cos_phi * y[3*idx] + cos_eps * sin_phi * y[3*idx+1] + sin_eps * y[3*idx + 2]);
            }
            else if ((T1 > 1e-5) && (T2 < 1e-5)) {
              dydx[3*idx] = 2.*(-xi * y[3*idx + 1] + D_Im * y[3*idx + 2]) - sin_eps * cos_phi / T1 * (sin_eps * cos_phi * y[3*idx] - sin_eps * sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq_final)));

              dydx[3*idx + 1] = 2.*(D_Re*y[3*idx + 2] + xi * y[3*idx]) + sin_phi * sin_eps / T1 * (sin_eps*cos_phi * y[3*idx]- sin_eps*sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq_final)));

              dydx[3*idx + 2] = -2.*(D_Im * y[3*idx] + D_Re * y[3*idx + 1]) - cos_eps / T1 * (sin_eps * cos_phi * y[3*idx]- sin_eps * sin_phi * y[3*idx+1] + cos_eps * y[3*idx + 2] - 1./2.*tanh(E_bog/(2. * T_eq_final)));
            }
            else {
                dydx[3*idx] = 2.*(-xi * y[3*idx+1] + D_Im * y[3*idx + 2]);

                dydx[3*idx + 1] = 2.*(D_Re * y[3*idx+2] + xi * y[3*idx]);

                dydx[3*idx + 2] = 2.*(-D_Im * y[3*idx] - D_Re * y[3*idx + 1]);
            }
        } // closes the double for loop
    } // closes void operator ()
};

int main(int argc, const char * argv[]) {
    cout << endl;
    clock_t t_total; // to measure runtime of simulation
    t_total = clock(); // start the clock t_total
    int L_kxy;
    double x_max;
     // Define the parameters:
    if (argc == 11){
      L_kxy = atof(argv[1]);
      L_kx = atof(argv[2]);
      N = L_kxy * L_kx;
      L_time = 200;
      t_min = 0.; // t_min of time evolution
      t_max = 20.; // t_max of time evolution
      J = 1.; // hopping parameter on square lattice
      mu = atof(argv[3]);//-1.18992; // chemical potential for filling = 0.5 and T=0
      V0 = 3; // strength of pairing interaction
      omega_D = J/2.;
      rho0 = 1/(2*pi*pi*J)*boost::math::ellint_1(sqrt(1-pow((mu)/(2*J), 2)/4));//0.135133;
      Delta0 = zbrent(gap_eq, 0.0001, 1.8, 1e-8);
      cout << "Delta0 = " << Delta0 << endl;
      temperature_init = 0.04 * Delta0;//0.453546 * 0.08547958631092525; // temperature of initial state
      A0 = sqrt(atof(argv[4])*Delta0);
      A0x = cos(0)*A0; // for momentum space calculation
      A0y = sin(0)*A0; // for momentum space calculation
      cout << "A0 = " << A0 << endl;
      tau = atof(argv[5])/Delta0;
      sigma = atof(argv[6])/Delta0;
      //sigma = tau/5;
      omega_pump = atof(argv[7]);//0.120865;
      T1 = atof(argv[8])*tau; //10./Delta0;
      T2 = atof(argv[9])*tau;
      x_max = atof(argv[10])*tau;
      cout << "Input parameters: L_kxy=" << L_kxy << ", L_kx=" << L_kx << ", A0=" << A0 << ", tau=" << tau << ", omega_pump=" << omega_pump << ", T1=" << T1 << ", T2=" << T2 << ", x_max=" << x_max << endl;
    }
    else{
      cout << "Default parameters" << endl << endl;
      L_kxy = 100;//number of discretization steps in momentum space
      L_kx = 100;
      N = L_kxy * L_kx;
      L_time = 200;

      t_min = 0.; // t_min of time evolution
      t_max = 20.; // t_max of time evolution
      J = 1.; // hopping parameter on square lattice
      mu = -1.18992; // chemical potential for filling = 0.5 and T=0
      V0 = 3; // strength of pairing interaction
      omega_D = J/2.;
      rho0 = 1/(2*pi*pi*J)*boost::math::ellint_1(sqrt(1-pow((mu)/(2*J), 2)/4));//0.135133;
      //Delta0 = 0.06076546880752605; // initial SC gap calculated from kxy-ky discretization
      //Delta0 = 0.08246266212661453; // initial SC gap calculated from kx-ky discretization
      Delta0 = zbrent(gap_eq, 0.0001, 1.8, 1e-8);
      cout << "Delta0 = " << Delta0 << endl;
      //Delta0 = 0.08547958631092525; //* 1.74 * sqrt(1 - 0.8); //0.0843798; // initial SC gap calculated from energy discretization at T=0 for V0 = 5, omega_D = J/2
      temperature_init = 0.04 * Delta0;//0.453546 * 0.08547958631092525; // temperature of initial state
      A0 = sqrt(1.5*Delta0);//0.410804;//sqrt(0.75*Delta0); //0.205946;//sqrt(0.04*Delta0); // amplitude of vector potential
      A0x = cos(0)*A0; // for momentum space calculation
      A0y = sin(0)*A0; // for momentum space calculation
      cout << "A0 = " << A0 << endl;
      delta = 0.; // phase shift between x and y components of laser field

      tau = 2*pi/Delta0;//59.2559;// 2 * pi * 5 / Delta0;//59.2559;//372.316; //* 0.120865 * 1.2 / Delta0;//5./Delta0; //2.*pi/Delta0; // width of pulse window
      sigma = tau/5;//pi/(2*Delta0);//tau/5.;//tau/5.; //pi/(4.*Delta0); // Gaussian width of pulse
      omega_pump = 0.120865;//0.120865;//2.5 * Delta0;//0.120865;//0.120865;;;//10./6.*Delta0;//0.759418;//0.120865;//0.759418;//0.120865;//10/6 * Delta0;//0.120865;//9.*Delta0; //3.*pi*Delta0; // center frequency of pump
      T1 = 0*tau; //10./Delta0;
      T2 = 0*tau; //10.*tau; //10./Delta0;
      x_max = 10*tau;
    }
    /*******************************************************************
     Define the momentum grid
     ********************************************************************/
    VecDoub kstates(3*L_kxy*L_kx);
    for(int idxy = 0; idxy < L_kxy; idxy++){
        double kxy = acos((omega_D-mu)/(2*J)-cos(A0))+(2*acos((omega_D+mu)/(-4*J*cos(A0))-0.5*(1/cos(A0)-1))-acos((omega_D-mu)/(2*J)-cos(A0)))/(static_cast<double>(L_kxy-1))*static_cast<double>(idxy);
        for (int idy = 0.; idy < L_kx; idy++){
            double ky = 2*acos((omega_D+mu)/(-4*J*cos(A0))-0.5*(1/cos(A0)-1))/(static_cast<double>(L_kx-1))*static_cast<double>(idy);
            if (kxy-ky>=0){
               kstates[3*idxy*L_kx+3*idy] = - 2*J * (cos(kxy-ky) + cos(ky)) - mu;
               kstates[3*idxy*L_kx+3*idy+1] = kxy - ky;
               kstates[3*idxy*L_kx+3*idy+2] = ky;
            }
            else {
               kstates[3*idxy*L_kx+3*idy] = -mu;
               kstates[3*idxy*L_kx+3*idy+1] = 0.;
               kstates[3*idxy*L_kx+3*idy+2] = 0.;
            }

        }
    }// end momentum grid
    double Delta0_k = gap_eq_sum_kxy(Delta0, L_kxy)+1;
    double ps_ratio = 1/Delta0_k;
    cout << "Delta0_k = " << Delta0_k << endl;
    cout << "phase space ratio = " << ps_ratio << endl;

    // Output options:
    bool output_delta = true;
    //cout << "K(0) = "  << boost::math::ellint_1(0) << endl; // test elliptical integral
    /*****************************************************************************
     ************** Solve ODE while laser is on 0 < t < tau **********************
     *****************************************************************************/
    const Int nvar = 3*N;
    const Doub atol = 1.0e-8;
    const Doub rtol = atol;
    const Doub h1 = 0.01;
    const Doub hmin = 0.0;

    const Doub x1 = 0.0; // initial time
    const Doub x1_1 = tau/2.;
    const Doub x1_2 = tau;
    const Doub x1_3 = 3.*tau/2.;
    const Doub x2 = 2.*tau; // final time

    VecDoub ystart(nvar);
    // Define the momentum grid: kstates
    for (int idx = 0; idx < N; idx++){
        double xi = kstates[3*idx];
        if (xi>=-omega_D && xi<=omega_D){
            ystart[3*idx] = 1./2.*Delta0/sqrt(pow(Delta0,2) + pow(xi,2)) * tanh(sqrt(pow(Delta0,2) + pow(xi,2))/(2*temperature_init));
            ystart[3*idx+1] = 0.;
            ystart[3*idx+2] = -1./2.*xi/sqrt(pow(Delta0,2) + pow(xi,2)) * tanh(sqrt(pow(Delta0,2)+pow(xi,2))/(2*temperature_init));
        }
        else {
            ystart[3*idx] = 0.;
            ystart[3*idx+1] = 0.;
            ystart[3*idx+2] = -1./2.*tanh(xi/(2*temperature_init));
        }
    }
    Output out1(500);
    Output out2(500);
    Output out3(500);
    Output out4(500);
    cout << "momentun space discretization used." << endl;
    rhs_bcs_k_space bcs_k(N, V0, J, mu, omega_D, Delta0, A0x, A0y, tau, sigma, omega_pump, T1, T2, kstates, ps_ratio);
    Odeint<StepperDopr853<rhs_bcs_k_space>> ode1(ystart, x1, x1_1, atol, rtol, h1, hmin, out1, bcs_k);
    ode1.integrate();

    for (int idx = 0; idx < N; idx++) {
        ystart[3*idx] = out1.ysave[3*idx][out1.count-1];
        ystart[3*idx+1] = out1.ysave[3*idx+1][out1.count-1];
        ystart[3*idx+2] = out1.ysave[3*idx+2][out1.count-1];
    }
    cout << "ode.integrate finished in region 0 < t < tau/2." << endl;
    Odeint<StepperDopr853<rhs_bcs_k_space>> ode2(ystart, x1_1, x1_2, atol, rtol, h1, hmin, out2, bcs_k);
    ode2.integrate();

    for (int idx = 0; idx < N; idx++) {
        ystart[3*idx] = out2.ysave[3*idx][out2.count-1];
        ystart[3*idx+1] = out2.ysave[3*idx+1][out2.count-1];
        ystart[3*idx+2] = out2.ysave[3*idx+2][out2.count-1];
    }
    cout << "ode.integrate finished in region tau/2 < t < tau." << endl;
    Odeint<StepperDopr853<rhs_bcs_k_space>> ode3(ystart, x1_2, x1_3, atol, rtol, h1, hmin, out3, bcs_k);
    ode3.integrate();

    for (int idx = 0; idx < N; idx++) {
        ystart[3*idx] = out3.ysave[3*idx][out3.count-1];
        ystart[3*idx+1] = out3.ysave[3*idx+1][out3.count-1];
        ystart[3*idx+2] = out3.ysave[3*idx+2][out3.count-1];
    }
    cout << "ode.integrate finished in region tau < t < 3*tau/2." << endl;
    Odeint<StepperDopr853<rhs_bcs_k_space>> ode4(ystart, x1_3, x2, atol, rtol, h1, hmin, out4, bcs_k);
    ode4.integrate();

    cout << "ode.integrate finished in region 3*tau/2 < t < 2tau. " << endl;
    // determine final energy in system when laser is off
    double DeltaRe_final = 0.;
    double DeltaIm_final = 0.;
    double energy_final = 0.;

    for (int jdx = 0; jdx < N; jdx++) {
        double xi = kstates[3*jdx];
        if (xi >= -omega_D && xi <= omega_D && kstates[3*jdx+1] >= 0){
            DeltaRe_final += ps_ratio * out4.ysave[3*jdx][out4.count-1]*(V0/N);
            DeltaIm_final += ps_ratio * out4.ysave[3*jdx+1][out4.count-1]*(-V0/N);
            energy_final += (ps_ratio/N)*(2.*xi*out4.ysave[3*jdx+2][out4.count-1] + sqrt(xi*xi + Delta0*Delta0));
        }
        //energy += (1/(N*N)) * (2.*xi*y[3*idx*N + 3*idy + 2] + sqrt(xi*xi + Delta0*Delta0));
    }
    energy_final += -1./V0*(DeltaRe_final*DeltaRe_final + DeltaIm_final*DeltaIm_final + Delta0*Delta0);
    energy_final = energy_final/(Delta0*Delta0);

    cout << "Energy = " << energy_final/rho0;

    double T_eq_final=0.;
    double Delta_eq_final=0.;

    double sol_u_final;
    double energy_final_TB = energy_final/rho0;

    if (energy_final_TB < 0.0001) {
        sol_u_final = 1/0.04;
        Delta_eq_final = 1. * Delta0;
        T_eq_final = 0.04 * Delta0;
    }
    else if (energy_final_TB > 1.557406319375) {
        sol_u_final = 0;
        Delta_eq_final = 0.;
        auto solve_0 = [&energy_final_TB](double T) {
            return solve_0_E(T, energy_final_TB);
        };
        T_eq_final = zbrent(solve_0, 0.56, 20, 1.0e-8) * Delta0; //0.566933
    }
    else {
        auto solve_u = [&energy_final_TB](double u) {
            return solve_u_E(u,energy_final_TB);
        };
        sol_u_final = zbrent(solve_u, 0.001, 5., 1.0e-8);
        Delta_eq_final = gap(sol_u_final) * Delta0;
        T_eq_final = Delta_eq_final / (2. * sol_u_final);
    }

    cout << "energy_final/Delta0^2 = " << energy_final << ", T_eq = " << T_eq_final << ", T/Tc = " << T_eq_final/(Delta0*0.566933) << ", Delta_eq/Delta0 = " << Delta_eq_final / Delta0 << endl;
    cout << " I am here" << endl;

    /**************************************************************************
     ************** Solve ODE while laser is off t > tau **********************
     **************************************************************************/
    const Doub x3 = 2.*tau; // initial time
    const Doub x3_1 = 5.*tau/2.;
    const Doub x4 = x_max; // final time

    for (int idx = 0; idx < N; idx++) {
        ystart[3*idx] = out4.ysave[3*idx][out4.count-1];
        ystart[3*idx+1] = out4.ysave[3*idx+1][out4.count-1];
        ystart[3*idx+2] = out4.ysave[3*idx+2][out4.count-1];
    }
    Output out_21(500);
    Output out_22(500);
    cout << "momentun space discretization used." << endl;

    rhs_bcs_k_space2 bcs_k_2(N, V0, J, mu, omega_D, Delta0, T_eq_final, Delta_eq_final, DeltaRe_final, DeltaIm_final, T1, T2, kstates, ps_ratio);
    Odeint<StepperDopr853<rhs_bcs_k_space2>> ode_21(ystart, x3, x3_1, atol, rtol, h1, hmin, out_21, bcs_k_2);
    ode_21.integrate();

    cout << " ode.integrate finished in region 2*tau < t < 5tau/2. " << endl;
    for (int idx = 0; idx < N; idx++) {
        ystart[3*idx] = out_21.ysave[3*idx][out_21.count-1];
        ystart[3*idx+1] = out_21.ysave[3*idx+1][out_21.count-1];
        ystart[3*idx+2] = out_21.ysave[3*idx+2][out_21.count-1];
    }
    Odeint<StepperDopr853<rhs_bcs_k_space2>> ode_22(ystart, x3_1, x4, atol, rtol, h1, hmin, out_22, bcs_k_2);
    ode_22.integrate();
    cout << " ode.integrate finished in region 5*tau/2 < t < 3tau. " << endl;
    /********************************************
     **************** Output ********************
     ********************************************/
    string L_kxy_string = to_string(L_kxy);
    string L_kx_string = to_string(L_kx);
    string A0x_string = to_string(A0x*A0x/Delta0);
    string A0y_string = to_string(A0y*A0y/Delta0);
    string sigma_string = to_string(sigma*Delta0);
    string tau_string = to_string(tau*Delta0);
    string omega_pump_string = to_string(omega_pump);
    string V0_string = to_string(V0);
    string T1_string = to_string(T1);
    string T2_string = to_string(T2);

    A0x_string.erase(A0x_string.find_last_not_of('0') + 1, string::npos);
    A0y_string.erase(A0y_string.find_last_not_of('0') + 1, string::npos);
    A0x_string.erase(A0x_string.find_last_not_of('.') + 1, string::npos);
    A0y_string.erase(A0y_string.find_last_not_of('.') + 1, string::npos);
    sigma_string.erase(sigma_string.find_last_not_of('0') + 1, string::npos);
    tau_string.erase(tau_string.find_last_not_of('0') + 1, string::npos);
    omega_pump_string.erase(omega_pump_string.find_last_not_of('0') + 1, string::npos);
    V0_string.erase(V0_string.find_last_not_of('0') + 1, string::npos);
    T1_string.erase(T1_string.find_last_not_of('0') + 1, string::npos);
    T2_string.erase(T2_string.find_last_not_of('0') + 1, string::npos);
    T1_string.erase(T1_string.find_last_not_of('.') + 1, string::npos);
    T2_string.erase(T2_string.find_last_not_of('.') + 1, string::npos);

    //outputFilename = "Output/Data-Momentum_states.dat";
    //outputFile.open(outputFilename);
    //outputFile.setf(ios::scientific);

    //for (int idx = 0; idx < kstates.size(); idx = idx + 3) {
    //    outputFile << kstates[idx] << " ";
    //    for (int jdx = 1; jdx <= 2; jdx++) {
    //        outputFile << kstates[idx + jdx] << " ";
    //    }
    //    outputFile << endl;
    //}
    //outputFile.close();


    if (output_delta == true) {
        outputFilename = "Output/Data(t)_len-Ax=" + A0x_string  + "-Ay=" + A0y_string + "-omega_p=" + omega_pump_string + "-tau=" + tau_string + "-sigma=" + sigma_string + "-T1=" + T1_string + "-T2=" + T2_string + "-L_kxy=" + L_kxy_string + "-L_kx=" + L_kx_string + ".dat";
        outputFile.open(outputFilename);

        outputFile.setf(ios::scientific);
        //double energy_out = 0.;
        //double energy_A_out = 0.;
        double Ax = 0.;
        double Ay = 0.;
        for (int idx = 0; idx < out1.count; idx++) {
            outputFile << out1.xsave[idx]*Delta0 << " ";
            double DeltaRe_out = 0.;
            double DeltaIm_out = 0.;
            double energy_out = 0.;
            if ((out1.xsave[idx] > 0.) && (out1.xsave[idx] < 2.*tau)) {
                Ax = A0x*exp(-pow(out1.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out1.xsave[idx]);
                Ay = A0y*exp(-pow(out1.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out1.xsave[idx]);
                //A = A0*exp(-pow(out.xsave[idx]-tau/2.,2)/(2.*sigma*sigma)) * cos(omega_pump*out.xsave[idx]);
            }
            for (int jdx = 0; jdx < N; jdx++) {
                //double xi = -2*J*(cos(kstates[3*jdx+1]) * cos(Ax) + cos(kstates[3*jdx+2])*cos(Ay))-mu;
                //if (xi>=-omega_D && xi<=omega_D && kstates[3*jdx+1] >= 0){

                    //energy_out += (ps_ratio/N) * (2.*xi*out.ysave[3*jdx+2][idx] + xi);
                //}
                if (kstates[3*jdx] >= -omega_D && kstates[3*jdx] <= omega_D && kstates[3*jdx+1]>=0){
                    DeltaRe_out += ps_ratio * out1.ysave[3*jdx][idx]*(V0/N);
                    DeltaIm_out += ps_ratio * out1.ysave[3*jdx+1][idx]*(-V0/N);
                    energy_out += (ps_ratio/N) * (2.*kstates[3*jdx]*out1.ysave[3*jdx+2][idx] + sqrt(kstates[3*jdx]*kstates[3*jdx] + Delta0*Delta0));
                }
                //if (kstates[3*jdx] >= -omega_D && kstates[3*jdx] <= omega_D && kstates[3*jdx+1]>=0){
                //    energy_out += (ps_ratio/N) * (2.*xi*out.ysave[3*jdx+2][idx] + xi - kstates[3*jdx]+ sqrt(kstates[3*jdx]*kstates[3*jdx] + Delta0*Delta0));
                //}
            }
            energy_out += -1./V0*(DeltaRe_out*DeltaRe_out + DeltaIm_out*DeltaIm_out + Delta0*Delta0);
            energy_out = energy_out/(Delta0*Delta0);
            outputFile << sqrt(DeltaRe_out*DeltaRe_out+DeltaIm_out*DeltaIm_out)/Delta0 << " " << energy_out << endl;
        }
        for (int idx = 0; idx < out2.count; idx++) {
            outputFile << out2.xsave[idx]*Delta0 << " ";
            double DeltaRe_out = 0.;
            double DeltaIm_out = 0.;
            double energy_out = 0.;
            if ((out2.xsave[idx] > 0.) && (out2.xsave[idx] < 2.*tau)) {
                Ax = A0x*exp(-pow(out2.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out2.xsave[idx]);
                Ay = A0y*exp(-pow(out2.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out2.xsave[idx]);
            }
            for (int jdx = 0; jdx < N; jdx++) {
                if (kstates[3*jdx] >= -omega_D && kstates[3*jdx] <= omega_D && kstates[3*jdx+1]>=0){
                    DeltaRe_out += ps_ratio * out2.ysave[3*jdx][idx]*(V0/N);
                    DeltaIm_out += ps_ratio * out2.ysave[3*jdx+1][idx]*(-V0/N);
                    energy_out += (ps_ratio/N) * (2.*kstates[3*jdx]*out2.ysave[3*jdx+2][idx] + sqrt(kstates[3*jdx]*kstates[3*jdx] + Delta0*Delta0));
                }
            }
            energy_out += -1./V0*(DeltaRe_out*DeltaRe_out + DeltaIm_out*DeltaIm_out + Delta0*Delta0);
            energy_out = energy_out/(Delta0*Delta0);
            outputFile << sqrt(DeltaRe_out*DeltaRe_out+DeltaIm_out*DeltaIm_out)/Delta0 << " " << energy_out << endl;
        }
        for (int idx = 0; idx < out3.count; idx++) {
            outputFile << out3.xsave[idx]*Delta0 << " ";
            double DeltaRe_out = 0.;
            double DeltaIm_out = 0.;
            double energy_out = 0.;
            if ((out3.xsave[idx] > 0.) && (out3.xsave[idx] < 2.*tau)) {
                Ax = A0x*exp(-pow(out3.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out3.xsave[idx]);
                Ay = A0y*exp(-pow(out3.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out3.xsave[idx]);
            }
            for (int jdx = 0; jdx < N; jdx++) {
                if (kstates[3*jdx] >= -omega_D && kstates[3*jdx] <= omega_D && kstates[3*jdx+1]>=0){
                    DeltaRe_out += ps_ratio * out3.ysave[3*jdx][idx]*(V0/N);
                    DeltaIm_out += ps_ratio * out3.ysave[3*jdx+1][idx]*(-V0/N);
                    energy_out += (ps_ratio/N) * (2.*kstates[3*jdx]*out3.ysave[3*jdx+2][idx] + sqrt(kstates[3*jdx]*kstates[3*jdx] + Delta0*Delta0));
                }
            }
            energy_out += -1./V0*(DeltaRe_out*DeltaRe_out + DeltaIm_out*DeltaIm_out + Delta0*Delta0);
            energy_out = energy_out/(Delta0*Delta0);
            outputFile << sqrt(DeltaRe_out*DeltaRe_out+DeltaIm_out*DeltaIm_out)/Delta0 << " " << energy_out << endl;
        }
        for (int idx = 0; idx < out4.count; idx++) {
            outputFile << out4.xsave[idx]*Delta0 << " ";
            double DeltaRe_out = 0.;
            double DeltaIm_out = 0.;
            double energy_out = 0.;
            if ((out4.xsave[idx] > 0.) && (out4.xsave[idx] < 2.*tau)) {
                Ax = A0x*exp(-pow(out4.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out4.xsave[idx]);
                Ay = A0y*exp(-pow(out4.xsave[idx]-tau,2)/(2.*sigma*sigma)) * cos(omega_pump*out4.xsave[idx]);
            }
            for (int jdx = 0; jdx < N; jdx++) {
                if (kstates[3*jdx] >= -omega_D && kstates[3*jdx] <= omega_D && kstates[3*jdx+1]>=0){
                    DeltaRe_out += ps_ratio * out4.ysave[3*jdx][idx]*(V0/N);
                    DeltaIm_out += ps_ratio * out4.ysave[3*jdx+1][idx]*(-V0/N);
                    energy_out += (ps_ratio/N) * (2.*kstates[3*jdx]*out4.ysave[3*jdx+2][idx] + sqrt(kstates[3*jdx]*kstates[3*jdx] + Delta0*Delta0));
                }
            }
            energy_out += -1./V0*(DeltaRe_out*DeltaRe_out + DeltaIm_out*DeltaIm_out + Delta0*Delta0);
            energy_out = energy_out/(Delta0*Delta0);
            outputFile << sqrt(DeltaRe_out*DeltaRe_out+DeltaIm_out*DeltaIm_out)/Delta0 << " " << energy_out << endl;
        }

        for (int idx = 0; idx < out_21.count; idx++) {
            outputFile << out_21.xsave[idx]*Delta0 << " ";
            double DeltaRe_out = 0.;
            double DeltaIm_out = 0.;
            double energy_out = 0.;
            for (int jdx = 0; jdx < N; jdx++) {
                double xi = kstates[3*jdx];
                if (xi >= -omega_D && xi <= omega_D && kstates[3*jdx+1] >= 0){
                    DeltaRe_out += ps_ratio * out_21.ysave[3*jdx][idx]*(V0/N);
                    DeltaIm_out += ps_ratio * out_21.ysave[3*jdx+1][idx]*(-V0/N);
                    energy_out += ps_ratio/N * (2.*xi*out_21.ysave[3*jdx+2][idx] + sqrt(xi*xi + Delta0*Delta0));
                }
            }
            energy_out += -1./V0*(DeltaRe_out*DeltaRe_out + DeltaIm_out*DeltaIm_out + Delta0*Delta0);
            energy_out = energy_out/(Delta0*Delta0);
            outputFile << sqrt(DeltaRe_out*DeltaRe_out+DeltaIm_out*DeltaIm_out)/Delta0 << " " << energy_out << endl;
        }
        for (int idx = 0; idx < out_22.count; idx++) {
            outputFile << out_22.xsave[idx]*Delta0 << " ";
            double DeltaRe_out = 0.;
            double DeltaIm_out = 0.;
            double energy_out = 0.;
            for (int jdx = 0; jdx < N; jdx++) {
                double xi = kstates[3*jdx];
                if (xi >= -omega_D && xi <= omega_D && kstates[3*jdx+1] >= 0){
                    DeltaRe_out += ps_ratio * out_22.ysave[3*jdx][idx]*(V0/N);
                    DeltaIm_out += ps_ratio * out_22.ysave[3*jdx+1][idx]*(-V0/N);
                    energy_out += ps_ratio/N * (2.*xi*out_22.ysave[3*jdx+2][idx] + sqrt(xi*xi + Delta0*Delta0));
                }
            }
            energy_out += -1./V0*(DeltaRe_out*DeltaRe_out + DeltaIm_out*DeltaIm_out + Delta0*Delta0);
            energy_out = energy_out/(Delta0*Delta0);
            outputFile << sqrt(DeltaRe_out*DeltaRe_out+DeltaIm_out*DeltaIm_out)/Delta0 << " " << energy_out << endl;
        }
        outputFile.close();
    }
    // Computational time:
    t_total = clock() - t_total; // time difference between current clock and clock at beginning of main
    cout << endl << "Actual total runtime = " << static_cast<double>(t_total)/CLOCKS_PER_SEC << " sec = " << static_cast<double>(t_total)/CLOCKS_PER_SEC/60. << " min = " << static_cast<double>(t_total)/CLOCKS_PER_SEC/60./60. << " h" << endl;
    cout << endl;
    return 0;
}
