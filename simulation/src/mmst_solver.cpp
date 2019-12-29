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



//#define ARMA_ALLOW_FAKE_GCC

#include <iostream>
#include <armadillo>
#include <cmath>
#include <complex>
#include <random>
#include <tuple>
#include "parameters.hpp"
#include <vector>
#include <string>

#ifdef __MPI_ENVIRONMENT__
#include <mpi.h>
#endif

#ifdef __OMP_ENVIRONMENT__
#include <omp.h>
#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#endif

using namespace std;






/*****************************************************************/
/* Definition of the 1D Cavity with Fundamental Functions        */
/*****************************************************************/

class OneDimCavity
{
private:
    size_t N, Ne, Nmodes;
    arma::cx_mat C, dCdt;
    arma::cx_mat Hs, HatP, SigmaX, SigmaZ;
    arma::cx_mat expHsdt, Identity22, expHsdtHalf;
    arma::mat OneMatrix21;
    arma::uvec PosiGridArry;
    double a, kFGR, kFGR_vacuum;
    arma::vec kFGR_array;
    arma::mat Xsi;

    size_t NGrids;
    double xmin, xmax, dx;
    arma::vec PosiCavity, PosiMolecules;
    arma::vec Ez, By, Jz, dEzdt, dBydt;

    double self_inteference_length, intdER2, intdBR2;
    arma::mat dER_mat, dBR_mat;
    
    // RK4 parameters
    arma::cx_mat dCdt1, dCdt2, dCdt3, dCdt4;
    arma::vec dEzdt1, dEzdt2, dEzdt3, dEzdt4;
    arma::vec dBydt1, dBydt2, dBydt3, dBydt4;
    // simulation and print parameters
    size_t ncount_C, ncount_EM, nprint_C, nprint_EM;
    double dt, tmax;

    arma::vec AlphaList, OmegaList, Pa, Qa;
public:
    // outputable parameters
    arma::mat traj_popu, traj_Ez, traj_Ez2; 
    
    OneDimCavity() : N(NMOLECULES), Ne(2), NGrids(NGRIDS), xmin(0.0), xmax(LCAVITY), ncount_C(0), ncount_EM(0), nprint_C(NPRINT_C), nprint_EM(NPRINT_EM), dt(DT), tmax(TMAX) 
    {
	/* Initialize parameters for wave functions */
	C = dCdt = arma::zeros<arma::cx_mat>(Ne, N);
	dCdt1 = dCdt2 = dCdt3 = dCdt4 = C;
	// set excited state distribution
	if (N < NMOLECULES_MAX)
	{
	    for (size_t i = 0; i < N; i++)
		C(1, i) = sqrt(EXCITEDSTATEDISTRIBUTION[i]);
	}
	else
	{
	    C.row(1).zeros();
	    C(1, int((N-1)/2)) = sqrt(EXCITEDSTATEDISTRIBUTION[0]);
	}
	//cout << "WaveFunction =\n" << arma::abs(C.row(1)) << endl;
	C.row(0) = 1.0 - C.row(1) % arma::conj(C.row(1));

	Hs = SigmaX = SigmaZ = arma::zeros<arma::cx_mat>(Ne, Ne);
	Hs(1, 1) = OMEGA0;

	SigmaX(0, 1) = SigmaX(1, 0) = 1.0;
	SigmaZ(0, 0) = 1.0;
	SigmaZ(1, 1) = -1.0;

	HatP = arma::zeros<arma::cx_mat>(Ne, Ne);
	HatP(0, 1) = HatP(1, 0) = 1.0;
	// split operator propagator exp(-iHdt) = exp(-iVdt/2)exp(-iHsdt)exp(-iVdt/2)
	expHsdt = arma::expmat(-I*Hs*dt);	
	expHsdtHalf = arma::expmat(-I*Hs*dt/2.0);	
	Identity22 = arma::eye<arma::cx_mat>(Ne, Ne);
	// the No. of grids of all TLSs
	PosiGridArry = arma::zeros<arma::uvec>(N);
	if (N < NMOLECULES_MAX)
	{
	    for(size_t i = 0; i < N; i++)
		PosiGridArry(i) = POSITIONGRIDARRAY[i];
	}
	else
	{
	    for(size_t i = 0; i < N; i++)
		PosiGridArry(i) = POSITIONGRIDARRAY[0] - ((N-1)/2 - i)*LATTICE_SPACING_GRIDS;
	}
	//cout << PosiGridArry << endl;

	kFGR = kFGR_vacuum = MU12 * MU12 * OMEGA0;
	kFGR_array = kFGR_vacuum * arma::ones<arma::vec>(N);	
	
	// This quantity will be used in Derivative() member function 
	OneMatrix21 = arma::ones<arma::mat>(Ne, 1);

	/* Initialize parameterts for EM fields in a cavity */
	PosiCavity = arma::linspace<arma::vec>(xmin, xmax, NGrids);
	
	Ez = By = Jz = dEzdt = dBydt = arma::zeros<arma::vec>(NGrids);
	dEzdt1 = dEzdt2 = dEzdt3 = dEzdt4 = Ez;
	dBydt1 = dBydt2 = dBydt3 = dBydt4 = By;

	dx = abs(PosiCavity(1) - PosiCavity(0));

	/* Parameter for the position of molecules */
	PosiMolecules = arma::zeros<arma::vec>(N);
	for(size_t i = 0; i < N; i++)
	    PosiMolecules(i) = PosiCavity(PosiGridArry(i));
	
	/* Initialize the polarization density distributions for molecules */
	a = 1.0 / 2.0 / pow(MOLECULARWIDTH, 2.0);
	Xsi = arma::zeros<arma::mat>(NGrids, N);
	for(size_t i = 0; i < N; i++)
	{
	    arma::vec X = PosiCavity - PosiMolecules(i);
	    Xsi.col(i) = ( MU12 * sqrt(a / PI) ) * arma::exp(-a * (X % X));
	}


	/* Parameters for output */
	traj_popu = arma::zeros<arma::mat>(nprint_C+1, Ne*N+1);
	traj_Ez   = traj_Ez2 = arma::zeros<arma::mat>(NGrids, nprint_EM+2);
	traj_Ez.col(0) = PosiCavity;	
	traj_Ez2.col(0) = PosiCavity;

	if (COHERENT_FIELD_AMPLITUDE > 0.0)
	    add_coherent_field();

	/* Parameters for photon modes */
	Nmodes = NMODES;
	AlphaList = arma::zeros<arma::vec>(Nmodes);
	for(size_t i = 0; i < Nmodes; i++)
	    AlphaList(i) = i+1;
	OmegaList = PI / LCAVITY * AlphaList;
	Qa = arma::zeros<arma::vec>(Nmodes);
	Pa = arma::zeros<arma::vec>(Nmodes);
    }
    
    void add_coherent_field()
    {
	double q_a = sqrt(1.0 / 2.0 / COHERENT_FIELD_OMEGA) * 2.0 * COHERENT_FIELD_AMPLITUDE * cos(PI / 8.0);
	arma::vec Ecoherent = sqrt(2.0 / LCAVITY) * COHERENT_FIELD_OMEGA * arma::sin(COHERENT_FIELD_OMEGA * PosiCavity) * q_a;
	Ez += Ecoherent;
	
	double p_a = -sqrt(COHERENT_FIELD_OMEGA / 2.0) * 2.0 * COHERENT_FIELD_AMPLITUDE * sin(PI / 8.0);
	p_a -= (DT / 2.0) * COHERENT_FIELD_OMEGA * COHERENT_FIELD_OMEGA * q_a;
	arma::vec Bcoherent = sqrt(2.0 / LCAVITY) * arma::cos(COHERENT_FIELD_OMEGA * (PosiCavity + dx / 2.0)) * p_a;
	By += Ecoherent;
    }

    void set_first_molecule_coherent_state(double p_ee = 0.1)
    {
	C(0, 0) = sqrt(1.0 - p_ee);
	C(1, 0) = sqrt(p_ee);
    }


    double get_random_01()
    {
	std::random_device rd; 
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> q(0.0, 1.0); 
	return q(gen);
    }
    
    double sgn(double ix)
    {
	return ( 0 < ix ) - (ix < 0);
    }

    // FDTD simulation on a Yee's cell
    void Step_FDTD(double t)
    {
	/* 1. calculate the derivative of wavefunction */
	// get the light-matter couplings at locations of TLSs
	for(size_t i = 0; i < N; i++)
	{   
	    double intEP = dx * arma::accu(Xsi.col(i) % Ez) ;
	    
	    double pre_sigmaz = OMEGA0 * dt / 2.0, pre_sigmax = intEP * dt;
	    double a = sqrt(pre_sigmax * pre_sigmax + pre_sigmaz * pre_sigmaz);
	    arma::cx_mat expHsc = arma::eye<arma::cx_mat>(Ne, Ne) * cos(a) + sin(a) * I * (pre_sigmax * SigmaX / a + pre_sigmaz * SigmaZ / a);
	    C.col(i) = expHsc * C.col(i);
	}

	/* 2. calculate current density according to the current wavefunction */
	Jz.zeros();
	arma::vec ImRho12 = arma::imag( C.row(0).st() % arma::conj(C.row(1).st()) );
	for(size_t i = 0; i < N; i++)
	    Jz += (-2.0 * OMEGA0 * (ImRho12(i)) ) * Xsi.col(i);

	/* 3. Propagate EM fields by FDTD algorithm */
	Ez.subvec(1, NGrids-1) += (dt / dx) * (By.subvec(0, NGrids-2) - By.subvec(1, NGrids-1));
	
	// add source
	Ez -= Jz * dt;
	
	// Boundary condition
	Ez(0) = Ez(NGrids-1) = 0.0;

	By.subvec(0, NGrids-2) += (dt / dx) * (Ez.subvec(0, NGrids-2) - Ez.subvec(1, NGrids-1));
    }

    // Propagate EM fields as photon modes
    void translate_modes_Ez()
    {
	for(size_t i = 0; i < NGrids; i++)
	    Ez(i) = sqrt(2.0 / LCAVITY) * arma::accu( OmegaList % arma::sin( PosiCavity(i)  * OmegaList) % Qa);
    }
    void Step_PhotonModes(double t)
    {
	/* 1. Calculate Light-Matter Interactions and propagate wave function for one time step*/
	arma::cx_mat C_old = C;
	for(size_t i = 0; i < N; i++)
	{  
	    double Ez_this_location = sqrt(2.0 / LCAVITY) * arma::accu( OmegaList % arma::sin( PosiMolecules(i)  * OmegaList) % Qa);
	    double intEP = MU12 * Ez_this_location; 
	    //arma::cx_mat expVdthalf = Identity22*cos(intEP*dt/2.0) + (I * sin(intEP*dt/2.0))*HatP;
	    //C.col(i) = expVdthalf * expHsdt * expVdthalf * C.col(i);

	    double pre_sigmaz = OMEGA0 * dt / 2.0, pre_sigmax = intEP * dt;
	    double a = sqrt(pre_sigmax * pre_sigmax + pre_sigmaz * pre_sigmaz);
	    arma::cx_mat expHsc = arma::eye<arma::cx_mat>(Ne, Ne) * cos(a) + (sin(a) * I) * (pre_sigmax * SigmaX / a + pre_sigmaz * SigmaZ / a);
	    C.col(i) = expHsc * C.col(i);
	}

	/* 2. Use split operator method to propagate photon modes for one time step*/
	arma::vec ReRho12 = arma::real( C_old.row(0).st() % arma::conj( C_old.row(1).st() ) );
	for (size_t j = ::START_FROM_JTH_MODE; j < Nmodes; j++)
	{
	    double dVdQj = 0.0;
	    for(size_t n = 0; n < N; n++)
		dVdQj -= MU12 * sqrt(2.0 / LCAVITY) * OmegaList(j) * sin( PosiMolecules(n)  * OmegaList(j) ) * (2.0 * ReRho12(n)) ;
	    Pa(j) -= dt / 2.0 * dVdQj;
	    double Pa_new = cos(OmegaList(j)*dt) * Pa(j) - OmegaList(j)*sin(OmegaList(j)*dt) * Qa(j); 
	    double Qa_new = 1.0 / OmegaList(j) * sin(OmegaList(j)*dt) * Pa(j) + cos(OmegaList(j)*dt) * Qa(j); 
	    Pa(j) = Pa_new;
	    Qa(j) = Qa_new;
	    Pa(j) -= dt / 2.0 * dVdQj;
	}

    }

    void Observe(double t, size_t icount, size_t nsteps_max)
    {
	// output traj
	if (icount % size_t(nsteps_max/nprint_C) == 0)
	{
	    traj_popu(ncount_C, 0) = t;
	    for(size_t i = 0; i < N; i++)
	    {
		traj_popu(ncount_C, 2*i+1) = pow(abs(C(0, i)), 2.0);
		traj_popu(ncount_C, 2*i+2) = pow(abs(C(1, i)), 2.0);
	    }
	    ncount_C++;
	}

	// output EM
	if (icount % size_t(nsteps_max/nprint_EM) == 0)
	{
	    //cout << "t = " << t << " No." << ncount_EM << " EM traj" << endl;
#ifdef RUN_MODES_INSTEAD_OF_FDTD
	    translate_modes_Ez();
#	endif
	    traj_Ez.col(ncount_EM+1) = Ez;
	    traj_Ez2.col(ncount_EM+1) = Ez % Ez;
	    ncount_EM++;
	}
    }

    friend class EhrenfestSingle;
};





/************************************************************************/
/* A Single Ehrenfest (Mean-Field) Simulation of 1D Cavity              */
/************************************************************************/

class EhrenfestSingle
{
private:
    bool  DoMultiTraj;
public:
    arma::vec E2_ZPE_correction, ZPE_Traj;
    OneDimCavity model;
    EhrenfestSingle():DoMultiTraj(false), model(OneDimCavity())
    {
	E2_ZPE_correction = arma::zeros<arma::vec>(model.NGrids);
    }

    void map_to_low_saturation(double max_population_scale=0.1)
    {
	cout << "Before scaling, excitated state population = \n" << arma::pow(arma::abs(model.C.row(1)), 2.0) << endl;
	model.C.row(1) *= max_population_scale;
	model.C.row(0) = arma::sqrt( 1.0 - arma::pow(model.C.row(1), 2.0) );
	cout << "After scaling,  excited state population = \n"   << arma::pow(arma::abs(model.C.row(1)), 2.0) << endl;
    }
    
    void sampling_wavefunction(double gamma = 0.42)
    {
	DoMultiTraj = true;
	for(size_t i = 0; i < model.N; i++)
	{
	    double p0 = pow( abs(model.C(0, i)), 2.0 );
	    double p1 = pow( abs(model.C(1, i)), 2.0 );

	    model.C(0, i) = sqrt(p0 + 2.0 * gamma * model.get_random_01() ) * exp(I*2.0*PI*model.get_random_01());
	    model.C(1, i) = sqrt(p1 + 2.0 * gamma * model.get_random_01() ) * exp(I*2.0*PI*model.get_random_01());
	    
	    //cout << "random number = " << model.get_random_01() << endl;
	}
    }
    

    void sampling_EM_vacuumFluctuations(size_t Nmodes)
    {
	DoMultiTraj = true;
	/* 1. Sampling over Q_a, P_a for a  = 0, 1, ..., Nmodes-1 */
	arma::vec AlphaList = arma::zeros<arma::vec>(Nmodes);
	for(size_t i = ::START_FROM_JTH_MODE; i < Nmodes; i++)
	    AlphaList(i) = i+1;
	arma::vec OmegaList = PI / LCAVITY * AlphaList;
	//cout << "OmegaList is \n" << OmegaList << endl;
	arma::vec Qa = arma::zeros<arma::vec>(Nmodes);
	arma::vec Pa = arma::zeros<arma::vec>(Nmodes);
	// assign data for p_a and q_a according to ZPE distribution of vacuum photonic states
	// the ZEP distribution should obey: rho(p_a, q_a) ~ 1/pi * exp(-p_a^2 / omega_a - omega_a * q_a^2)
	std::random_device rd;
	std::mt19937 gen(rd());
	for(size_t i = ::START_FROM_JTH_MODE; i < Nmodes; i++)
	{
	    double variance_p = sqrt( OmegaList(i) / 2.0 );
	    std::normal_distribution<double> p_alpha(0.0, variance_p);
	    Pa(i) = p_alpha(gen);

	    double variance_q = sqrt( 1.0 / OmegaList(i) / 2.0 );
	    std::normal_distribution<double> q_alpha(0.0, variance_q);
	    Qa(i) = q_alpha(gen);
	}

	model.Pa = Pa;
	model.Qa = Qa;

	arma::vec Pa_half = Pa - (model.dt / 2.0) * OmegaList % OmegaList % Qa;
	/* 2. Transform (P, Q) to (E, B) according to the following transformation */
	// E(r, t) = sum_a sqrt(2 / Lcavity) * omega_a * sin(a*pi*r/Lcavity)
	for(size_t i = 0; i < model.NGrids; i++)
	{
	    model.Ez(i) += sqrt(2.0 / LCAVITY) * arma::accu( OmegaList % arma::sin( model.PosiCavity(i)  * OmegaList) % Qa);
	    //model.By(i) += sqrt(2.0 / LCAVITY) * arma::accu(             arma::cos( model.PosiCavity(i)  * OmegaList) % Pa);
	    // I need to consider change the initialization of By, becasuse I propagate E,B in the Yee's cell, so they should
	    // have some little difference instead of the expression in the same (t, x). A slight change may lead to some big
	    // change. This needs to be investigated!
	    model.By(i) += sqrt(2.0 / LCAVITY) * arma::accu(             arma::cos( (model.PosiCavity(i) + model.dx/2.0)  * OmegaList) % Pa_half);
	}

	/* 3. Calculate the E2_ZPE_correction if necessary */
	E2_ZPE_correction = arma::zeros<arma::vec>(model.NGrids);
	for(size_t i = 0; i < model.NGrids; i++)
	{
	    E2_ZPE_correction(i) = arma::sum( arma::pow(arma::sqrt(OmegaList / LCAVITY) % arma::sin( model.PosiCavity(i) * OmegaList ), 2.0) );
	}
    }
    

    void sampling_EM_vacuumFluctuations_circle(size_t Nmodes)
    {
	DoMultiTraj = true;
	/* 1. Sampling over Q_a, P_a for a  = 0, 1, ..., Nmodes-1 */
	arma::vec AlphaList = arma::zeros<arma::vec>(Nmodes);
	for(size_t i = 0; i < Nmodes; i++)
	    AlphaList(i) = i+1;
	arma::vec OmegaList = PI / LCAVITY * AlphaList;
	//cout << "OmegaList is \n" << OmegaList << endl;
	arma::vec Qa = arma::zeros<arma::vec>(Nmodes);
	arma::vec Pa = arma::zeros<arma::vec>(Nmodes);
	// assign data for p_a and q_a according to ZPE distribution of vacuum photonic states
	// the ZEP distribution should obey: rho(p_a, q_a) ~ 1/pi * exp(-p_a^2 / omega_a - omega_a * q_a^2)
	std::random_device rd;
	std::mt19937 gen(rd());
	ZPE_Traj = arma::zeros<arma::vec>(Nmodes);

	for(size_t i = 0; i < Nmodes; i++)
	{
	    std::normal_distribution<double> random_gaussian(0.0, 1.0);
	    double temp = random_gaussian(gen);
	    double random_angle = 2.0 * PI * model.get_random_01();

	    temp = 1.0 / 8.0;
	    double variance_p = sqrt( OmegaList(i) );
	    Pa(i) = variance_p * temp * cos(random_angle);
	    double variance_q = sqrt( 1.0 / OmegaList(i) );
	    Qa(i) = variance_q * temp * sin(random_angle);

	    ZPE_Traj(i) = pow(Pa(i), 2.0) + pow(OmegaList(i) * Qa(i), 2.0);
	}

	model.Pa = Pa;
	model.Qa = Qa;

	arma::vec Pa_half = Pa - (model.dt / 2.0) * OmegaList % OmegaList % Qa;
	/* 2. Transform (P, Q) to (E, B) according to the following transformation */
	// E(r, t) = sum_a sqrt(2 / Lcavity) * omega_a * sin(a*pi*r/Lcavity)
	for(size_t i = 0; i < model.NGrids; i++)
	{
	    model.Ez(i) += sqrt(2.0 / LCAVITY) * arma::accu( OmegaList % arma::sin( model.PosiCavity(i)  * OmegaList) % Qa);
	    model.By(i) += sqrt(2.0 / LCAVITY) * arma::accu(             arma::cos( (model.PosiCavity(i) + model.dx/2.0)  * OmegaList) % Pa_half);
	}

	/* 3. Calculate the E2_ZPE_correction if necessary */
	E2_ZPE_correction = arma::zeros<arma::vec>(model.NGrids);
	for(size_t i = 0; i < model.NGrids; i++)
	{
	    E2_ZPE_correction(i) = arma::sum( arma::pow(arma::sqrt(OmegaList / LCAVITY) % arma::sin( model.PosiCavity(i) * OmegaList ), 2.0) );
	}
    }
    

    void change_sign_circle(size_t Nmodes)
    {
	for(size_t i = 0; i < Nmodes; i++)
	{
	    double length = pow( model.Pa(i), 2.0) + pow(model.OmegaList(i) * model.Qa(i), 2.0);
	    if (length < ZPE_Traj(i)*0.0)
		model.Pa(i) *= -1.0;
	}
    }

    std::tuple<arma::mat, arma::mat, arma::mat> run()
    {
	// Do simulation
	double dt = model.dt, tmax = model.tmax;
	size_t ncount = 0, nstepmax = size_t(tmax/dt+1);
	for(double t = 0.0; t < tmax; t+=dt)
	{
	    model.Observe(t, ncount, nstepmax);
	    model.Step_FDTD(t);
	    ncount++;
	}
	// Finally, we output results
	if (!DoMultiTraj)
	{
	    model.traj_popu.save("traj_Eh.txt", arma::raw_ascii);
	    model.traj_Ez.save("Ez_Eh.txt", arma::raw_ascii);
	    model.traj_Ez2.save("Ez2_Eh.txt", arma::raw_ascii);
	}

	return std::make_tuple(model.traj_popu, model.traj_Ez, model.traj_Ez2);
    }
    
    std::tuple<arma::mat, arma::mat, arma::mat> run_modes()
    {
	// Do simulation
	double dt = model.dt, tmax = model.tmax;
	size_t ncount = 0, nstepmax = size_t(tmax/dt+1);
	for(double t = 0.0; t < tmax; t+=dt)
	{
	    model.Observe(t, ncount, nstepmax);
	    model.Step_PhotonModes(t);
	    //change_sign_circle(NMODES);
	    ncount++;
	}
	// Finally, we output results
	if (!DoMultiTraj)
	{
	    model.traj_popu.save("traj_Eh.txt", arma::raw_ascii);
	    model.traj_Ez.save("Ez_Eh.txt", arma::raw_ascii);
	    model.traj_Ez2.save("Ez2_Eh.txt", arma::raw_ascii);
	}

	return std::make_tuple(model.traj_popu, model.traj_Ez, model.traj_Ez2);
    }
};




/************************************************************************/
/* MultiTrajectory Ehrenfest Dynamics of 1D Cavity                      */
/************************************************************************/

class EhrenfestMulti
{
private:
    size_t ntraj;
    double gamma;
    bool sampling_C, sampling_EM;

    arma::mat RecordWaveFunctionDistr;
    bool record_C;
    size_t NInstanceMax;

    std::string method_name;
public:
    EhrenfestMulti(size_t intraj=1280, bool isampling_C=true, bool isampling_EM=true, bool irecord_C=false):ntraj(intraj), gamma(GAMMA_ELECTRONIC_ZPE),
		    sampling_C(isampling_C), sampling_EM(isampling_EM), record_C(irecord_C)
    {
	if (NTRAJ_MULTITRAJ_Eh > 512)
	    NInstanceMax  = 512;
	else NInstanceMax = NTRAJ_MULTITRAJ_Eh;

	// construct a recorder to record the wave function distribution at disered place, save up to NInstanceMax instances
	if (record_C)
	    RecordWaveFunctionDistr = arma::zeros<arma::mat>(NInstanceMax, 2*NMOLECULES*(NPRINT_EM+1));   
	
	if (sampling_C && sampling_EM)
	    method_name = "MultiEh";
	else if (sampling_C && !sampling_EM)
	    method_name = "PreBinSQC";
	else if (!sampling_C && sampling_EM)
	    method_name = "StochasticED";
    }
    
    std::tuple<arma::mat, arma::mat, arma::mat> run_serial(int num_procs = 1, int id = 0)
    {
	// construct the averaged data
	arma::mat traj_popu_avg = arma::zeros<arma::mat>(NPRINT_C+1, 2*NMOLECULES+1);
	arma::mat traj_Ez_avg, traj_Ez2_avg;
	traj_Ez_avg = traj_Ez2_avg = arma::zeros<arma::mat>(NGRIDS, NPRINT_EM+2);

	size_t ntraj_thisProc = 0;	
	
#ifdef __OMP_ENVIRONMENT__
	#pragma omp parallel for shared(RecordWaveFunctionDistr, ntraj_thisProc) reduction(+ : traj_popu_avg, traj_Ez_avg, traj_Ez2_avg)	
#endif
	for(size_t n = id; n < ntraj; n+=num_procs)
	{
	    
	    EhrenfestSingle Eh;
	    // do sampling over zero-point energy of classical photon modes
	    if (sampling_EM) 
		Eh.sampling_EM_vacuumFluctuations(NMODES);
	    // do sampling over zero-point energy of electronic DOFs
	    if (sampling_C)
	    	Eh.sampling_wavefunction(gamma);
	    // run dynamics
#ifdef RUN_MODES_INSTEAD_OF_FDTD
	    auto results = Eh.run_modes();
#else
	    auto results = Eh.run();
#endif
	    cout << "END of No."<< n << " " << method_name << " traj" << endl;
	    
	    traj_popu_avg += std::get<0>(results);
	    traj_Ez_avg += std::get<1>(results);
	    arma::mat E2_raw = std::get<2>(results);
	    for(size_t i = 0; i < NPRINT_EM+1; i++)
	    	E2_raw.col(i+1) -= Eh.E2_ZPE_correction;
	    traj_Ez2_avg += E2_raw;

	    // finally, according to the traj of electronic DoFs, record distribution
	    if (record_C && n < NInstanceMax)
	    {
		size_t every_slice = size_t(double(NPRINT_C+1) / double(NPRINT_EM+1));
		arma::mat traj_popu = std::get<0>(results);
		for(size_t i = 0; i < NPRINT_EM+1; i++)
		{
		    arma::vec slice = traj_popu.row(i*every_slice).t();
		    if (sampling_C) 
			slice -= GAMMA_ELECTRONIC_ZPE;
		    RecordWaveFunctionDistr(n, arma::span(i*2*NMOLECULES, (i+1)*2*NMOLECULES-1)) = slice.subvec(1, 2*NMOLECULES).t();
		}
	    }

	    ntraj_thisProc +=1;
	}
	
	traj_popu_avg /= ntraj_thisProc;
	traj_Ez_avg   /= ntraj_thisProc;
	traj_Ez2_avg  /= ntraj_thisProc;

#ifndef __MPI_ENVIRONMENT__
	if (num_procs == 1 && id == 0)
	    save_data(traj_popu_avg, traj_Ez_avg, traj_Ez2_avg, RecordWaveFunctionDistr);
#endif	
	return std::make_tuple(traj_popu_avg, traj_Ez_avg, traj_Ez2_avg);
    }

    void run_parallel(int num_procs, int id)
    {
	// Get the averaged trajs for each processor
	auto results = run_serial(num_procs, id);
	arma::mat traj_popu_avg = std::get<0>(results);
	arma::mat traj_Ez_avg = std::get<1>(results);
	arma::mat traj_Ez2_avg = std::get<2>(results);

	// After calculation, we reduce these data together
	arma::mat traj_popu_avg_global, traj_Ez_avg_global, traj_Ez2_avg_global, RecordWaveFunctionDistr_global;
	
	if (id == 0)
	{
	    traj_popu_avg_global = arma::zeros<arma::mat>(arma::size(traj_popu_avg));
	    traj_Ez_avg_global = arma::zeros<arma::mat>(arma::size(traj_Ez_avg)); 
	    traj_Ez2_avg_global = arma::zeros<arma::mat>(arma::size(traj_Ez2_avg)); 
	    RecordWaveFunctionDistr_global = arma::zeros<arma::mat>(arma::size(RecordWaveFunctionDistr)); 
	}

#ifdef __MPI_ENVIRONMENT__
	::MPI_Barrier(MPI_COMM_WORLD);
	::MPI_Reduce(traj_popu_avg.memptr(), traj_popu_avg_global.memptr(), traj_popu_avg.n_elem, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	::MPI_Reduce(traj_Ez_avg.memptr(),   traj_Ez_avg_global.memptr(),   traj_Ez_avg.n_elem,   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	::MPI_Reduce(traj_Ez2_avg.memptr(),  traj_Ez2_avg_global.memptr(),  traj_Ez2_avg.n_elem,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	::MPI_Reduce(RecordWaveFunctionDistr.memptr(),  RecordWaveFunctionDistr_global.memptr(),  RecordWaveFunctionDistr.n_elem,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	::MPI_Barrier(MPI_COMM_WORLD);
#else
	traj_popu_avg_global = traj_popu_avg;
	traj_Ez_avg_global = traj_Ez_avg;
	traj_Ez2_avg_global = traj_Ez2_avg;
#endif
	traj_popu_avg_global /= num_procs;
	traj_Ez_avg_global /= num_procs;
	traj_Ez2_avg_global /= num_procs;

	if (id == 0)
	    save_data(traj_popu_avg_global, traj_Ez_avg_global, traj_Ez2_avg_global, RecordWaveFunctionDistr_global);
    }
    
    void save_data(arma::mat traj_popu_avg, arma::mat traj_Ez_avg, arma::mat traj_Ez2_avg, arma::mat RecordWaveFunctionDistr_tot)
    {
	if (sampling_C)
	    for(size_t i = 0; i < 2*NMOLECULES; i++)
		traj_popu_avg.col(i+1) -= GAMMA_ELECTRONIC_ZPE;

	// save data	
	traj_popu_avg.save("traj_" + method_name + ".txt", arma::raw_ascii);
	traj_Ez_avg.save("Ez_" + method_name + ".txt", arma::raw_ascii);
	traj_Ez2_avg.save("Ez2_" + method_name + ".txt", arma::raw_ascii);
	
	if (record_C)
	    RecordWaveFunctionDistr_tot.save("ElectronicDist_" + method_name + ".txt", arma::raw_ascii);
    }
};



int main()
{
    int num_procs = 1, id = 0;
    
#ifdef __MPI_ENVIRONMENT__
    ::MPI_Init(nullptr, nullptr);
    ::MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    ::MPI_Comm_rank(MPI_COMM_WORLD, &id);
#endif

    /* Simulation of multitraj Ehrenfest dynamics with parallel computing available*/ 
    EhrenfestMulti model_MT(NTRAJ_MULTITRAJ_Eh, ::SAMPLING_C, ::SAMPLING_EM, ::OUTPUT_ELECTRONIC_DIST_MTEh);
    model_MT.run_parallel(num_procs, id);
    
    //EhrenfestMulti model_SED(NTRAJ_MULTITRAJ_Eh, false, ::SAMPLING_EM, ::OUTPUT_ELECTRONIC_DIST_MTEh);
    //model_SED.run_parallel(num_procs, id);
    
    //EhrenfestMulti model_SQC(NTRAJ_MULTITRAJ_Eh, ::SAMPLING_C, false, ::OUTPUT_ELECTRONIC_DIST_MTEh);
    //model_SQC.run_parallel(num_procs, id);

#ifdef __MPI_ENVIRONMENT__    
    ::MPI_Finalize();
#endif

    return 0;
}
