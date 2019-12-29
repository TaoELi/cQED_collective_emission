/*****************************************************************/
/* Created and uploaded  by Tao E. Li @ 2019.12.28       
 * 
 * If this code is helpful for your research, please cite:
 *
 * Li, T. E., Chen, H.-T., Nitzan, A., & Subotnik, J. E. (2019). 
 * Quasiclassical Modeling of Cavity Quantum Electrodynamics. 
 * http://arxiv.org/abs/1910.02299 
 *
 * For any questions, please contact
 *  taoli@sas.upenn.edu
 * or
 *  t.e.li@outlook.com                                           */
/*****************************************************************/

#ifndef __PARAMETERS_HPP__
#define __PARAMETERS_HPP__

#define __MPI_ENVIRONMENT__
//#define __OMP_ENVIRONMENT__

// #define RUN_MODES_INSTEAD_OF_FDTD


#include <iostream>
#include <complex>
#include <armadillo>
#include <cmath>

using namespace std;

// Global variables
const double PI = arma::datum::pi;
const std::complex<double> I = std::complex<double>(0.0, 1.0);

const size_t NGRIDS = 5001;

const double LCAVITY = 2.0 * PI;

const size_t NMOLECULES = 1;
const size_t NMOLECULES_MAX = 3;
//const size_t POSITIONGRIDARRAY[] =        {2440, 2460, 2480, 2500, 2520, 2540, 2560};
//const double EXCITEDSTATEDISTRIBUTION[] = {0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0};
const size_t LATTICE_SPACING_GRIDS = 25;
const size_t POSITIONGRIDARRAY[] = {2500};
const double EXCITEDSTATEDISTRIBUTION[] = {1.0};


const double MU12 = sqrt(PI) * 0.1;
const double OMEGA0 = 100.0;
const double MOLECULARWIDTH = 1e-3;

const size_t NPRINT_C = 1000;
const size_t NPRINT_EM = 10;

const double DT = PI / (NGRIDS - 1);
const double TMAX = PI * 3.0;

const size_t NMODES = 400;

// Parameters for the added coherent E-field
const double COHERENT_FIELD_OMEGA = OMEGA0;
const double COHERENT_FIELD_AMPLITUDE = 0.0;

// parameters for Ehrenfest+R dynamics
const size_t NTRAJ_EHRENFEST_R = 128;
const bool OUTPUT_ELECTRONIC_DIST_EHR = true;

// parameters for multitrajectory Ehrenfest dynamics
const size_t NTRAJ_MULTITRAJ_Eh = 128;
const double GAMMA_ELECTRONIC_ZPE = 0.45;
const bool SAMPLING_C = true;
const bool SAMPLING_EM = true;
const bool OUTPUT_ELECTRONIC_DIST_MTEh = true;

const size_t START_FROM_JTH_MODE = 0;

#endif
