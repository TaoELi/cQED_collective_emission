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

using namespace std;


/************************************************************************/
/* Quantum Mechanical Simulation of 1D Cavity winthin CIS Approximation */
/************************************************************************/

class QuantumSimulator
{
private:
    // necessary parameters
    size_t Nmodes, Ndim, NGrids;
    arma::vec AlphaList, OmegaList;
    
    // wave function and Hamiltonians
    arma::cx_vec C_CIS;
    arma::cx_mat H_atoms, H_EM, H_int, H_CIS, expHCISdt;
    
    // Observe operators
    std::vector<arma::cx_mat> PExcitedS;
    std::vector<arma::cx_mat> EOperator, E2Operator;

    // output variables
    arma::mat traj_Ez, traj_Ez2, traj_popu;

    // simulation control parameters
    double dt, tmax;
    size_t ncount, nprint_C, nprint_EM, ncount_C, ncount_EM;

    // RK4 parameters
    arma::cx_vec dCdt1, dCdt2, dCdt3, dCdt4;
public:

    QuantumSimulator():Nmodes(NMODES), NGrids(400), dt(DT), tmax(TMAX), ncount(0), nprint_C(NPRINT_C), nprint_EM(NPRINT_EM), 
		       ncount_C(0), ncount_EM(0)
    {
	// 1. Construct the CIS wave function
	// |Pis(t)> = c0(t)|g>|g>_j|0>_k + c1(t) |e>|g>_j|0>_k + \sum_{j=2}^{M} cj(t) |g>|e>_j|0>k + \sum_{k} d_k(t)|g>|g>_j|1k>_k
	Ndim = NMOLECULES + NMODES + 1;

	size_t NEM_start = NMOLECULES + 1;
	size_t NMole_start = 1;

	C_CIS = arma::zeros<arma::cx_vec>(Ndim);
	// set the wave function; note that we only set the first molecule excited and ALL OTHER molecules are in ground state
	if (NMOLECULES < NMOLECULES_MAX)
	    for(size_t i = 0; i < NMOLECULES; i++)
	    {
		double pe = EXCITEDSTATEDISTRIBUTION[i];
		C_CIS(NMole_start+i) = sqrt(pe);
	    }
	else
	{
	    C_CIS.subvec(NMole_start, NMole_start+NMOLECULES-1).zeros();
	    C_CIS(NMole_start+size_t((NMOLECULES-1)/2)) = sqrt(EXCITEDSTATEDISTRIBUTION[0]);
	}
	C_CIS(0) = sqrt(1.0 - arma::accu(arma::pow(arma::abs(C_CIS.subvec(NMole_start, NMole_start+NMOLECULES-1)),2.0)) );
	//cout << C_CIS.subvec(NMole_start, NMole_start+NMOLECULES-1) << endl;	
	// 2. Construct the CIS Hamiltonian, H_CIS = H_Atoms + H_EM + H_int
	// 2.0 calculate the energies of photonic modes
	AlphaList = arma::zeros<arma::vec>(Nmodes);
	for(size_t i = 0; i < Nmodes; i++)
	    AlphaList(i) = i+1;
	OmegaList = PI / LCAVITY * AlphaList;

	// 2.1 H_EM: <0k | H_EM |0k' or 1k'> = 0.0;  <1k | H_EM |1k'> = omega_k * delta_kk'
	H_EM = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	for(size_t i = 0; i < Nmodes; i++)
	    H_EM(NEM_start+i, NEM_start+i) = OmegaList(i);
	
	// 2.2 H_atoms: add all molecules' omega_0 / 2.0 * SigmaZ, where SigmaZ = |e><e| - |g><g|
	H_atoms = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	for(size_t i = 0; i < NMOLECULES; i++)
	{
	    // build SigmaX for No.i molecule
	    arma::cx_mat SigmaZ = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	    SigmaZ.diag() += -1.0; // corresponding to -|g><g|
	    SigmaZ(NMole_start + i, NMole_start + i) = 1.0; // corresponding to |e><e|
	    H_atoms += (OMEGA0 / 2.0) * SigmaZ;
	}

	// 2.3 H_int: sum over both molecules and photon modes
	H_int = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	// positions of atoms
	arma::uvec PosiGridArry = arma::zeros<arma::uvec>(NMOLECULES);
	if (NMOLECULES < NMOLECULES_MAX)
	{
	    for(size_t i = 0; i < NMOLECULES; i++)
		PosiGridArry(i) = POSITIONGRIDARRAY[i];
	}
	else
	{
	    for(size_t i = 0; i < NMOLECULES; i++)
		PosiGridArry(i) = POSITIONGRIDARRAY[0] - ((NMOLECULES-1)/2 - i)*LATTICE_SPACING_GRIDS;
	}
	/* Initialize parameterts for EM fields in a cavity */
	arma::vec PosiCavity = arma::linspace<arma::vec>(0.0, LCAVITY, ::NGRIDS);
	/* Parameter for the position of molecules */
	arma::vec PosiMolecules = arma::zeros<arma::vec>(NMOLECULES);
	for(size_t i = 0; i < NMOLECULES; i++)
	    PosiMolecules(i) = PosiCavity(PosiGridArry(i));
	
	arma::vec AtomPositions = PosiMolecules;
	
	// create H_int finally       
	for(size_t i = 0; i < Nmodes; i++)
	{
	    for(size_t j = 0; j < NMOLECULES; j++)
	    {
		double gn_j = sqrt(OmegaList(i) / LCAVITY) * MU12 * sin(OmegaList(i) * AtomPositions(j));
		H_int(NMole_start+j, NEM_start+i) -= gn_j;
		H_int(NEM_start+i, NMole_start+j) -= gn_j;
	    }
	}

       // 2.4 Create the total Hamiltonian
       H_CIS = H_EM + H_atoms + H_int;
       expHCISdt = arma::expmat(-I*H_CIS*dt);       
       // 3. Finally, construct some operators for obervation
       for(size_t i = 0; i < NMOLECULES; i++)
       {
	   arma::cx_mat PExcited = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	   PExcited(NMole_start+i, NMole_start+i) = 1.0; // corresponding to |e><e|
	   PExcitedS.push_back(PExcited);
       }
       
       arma::vec X = arma::linspace<arma::vec>(0.0, LCAVITY, NGrids);	
       
       for(size_t i = 0; i < NGrids; i++)
       {
	   // prepare the E and E^2 operator for position X(i)
	   arma::cx_mat E = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	   arma::cx_mat E2 = arma::zeros<arma::cx_mat>(Ndim, Ndim);
	   
	   // we need to sum over modes to do summation for E and E^2
	   for(size_t j = 0; j < Nmodes; j++)
	   {
	       double coeff_j = sqrt( OmegaList(j) / LCAVITY )  * sin(OmegaList(j)*X(i));
	       E(0, NEM_start+j) += coeff_j;
	       E(NEM_start+j, 0) += coeff_j;
	       
	       // E^2 Operator = sum_j,k coeff_j * coeff_k [a_j + aDagger_j][a_k + aDagger_k]
	       // note that a|n> = \sqrt(n)|n-1>; aDagger|n> = \sqrt(n+1)|n+1>
	       // a_j * a_k, aDagger_j * aDagger_k are null in CIS basis; 
	       // aDagger_j * a_k has value
	       // Normal ordered a_j * aDagger_k = aDagger_k * a_j = aDagger_j * a_k
	       for(size_t k = 0; k < Nmodes; k++)
	       {
		   double coeff_k = sqrt( OmegaList(k) / LCAVITY ) * sin(OmegaList(k)*X(i));
		   // 3. Normal-ordered cross term = 2.0 * (aDagger_j * a_k)
		   // <1j|aDagger_j * a_k|1k> = 1.0
		   E2(NEM_start+j, NEM_start+k) += 2.0 * coeff_j * coeff_k;
	       }	   
	   }
	   // append operators to std::vector 
	   EOperator.push_back(E);
	   E2Operator.push_back(E2);
       }
       
       
       traj_Ez   = traj_Ez2 = arma::zeros<arma::mat>(NGrids, nprint_EM+2);
       traj_popu = arma::zeros<arma::mat>(NPRINT_C+1, 2*NMOLECULES+1);
       traj_Ez.col(0) = X;	
       traj_Ez2.col(0) = X;	
       cout << "End of initializer" << endl;

       // 4. Initialize RK4 parameters
       dCdt1 = dCdt2 = dCdt3 = dCdt4 = arma::zeros<arma::cx_vec>(Ndim);
    }

    void Observe(double t, size_t icount, size_t nsteps_max)
    {
	// output traj
	if (icount % size_t(nsteps_max/nprint_C) == 0)
	{
	    traj_popu(ncount_C, 0) = t;
	    for(size_t i = 0; i < NMOLECULES; i++)
	    {
		double pe = abs( arma::as_scalar(C_CIS.t() * PExcitedS[i] * C_CIS) );
		traj_popu(ncount_C, 2*i+1) = 1.0 - pe;
		traj_popu(ncount_C, 2*i+2) = pe;
	    }
	    ncount_C++;
	}

	// output EM
	if (icount % size_t(nsteps_max/nprint_EM) == 0)
	{
	    cout << "Outputing No." << ncount_EM << " EM" << endl;
	    arma::vec Ez = arma::zeros<arma::vec>(NGrids);
	    for(size_t i = 0; i < NGrids; i++)
		Ez(i) = std::real( arma::as_scalar(C_CIS.t() * EOperator[i] * C_CIS) );

	    traj_Ez.col(ncount_EM+1) = Ez;
	    
	    arma::vec Ez2 = arma::zeros<arma::vec>(NGrids);
	    for(size_t i = 0; i < NGrids; i++)
		Ez2(i) = std::abs( arma::as_scalar(C_CIS.t() * E2Operator[i] * C_CIS) );
	    traj_Ez2.col(ncount_EM+1) = Ez2;
	    ncount_EM++;
	}
    }

    void derivative(const arma::cx_vec &iC, arma::cx_vec &idCdt)
    {
	idCdt = -I * H_CIS * iC;
    }

    void Step_RK4()
    {
	derivative(C_CIS, dCdt1);
	derivative(C_CIS + (dt / 2.0) * dCdt1, dCdt2);
	derivative(C_CIS + (dt / 2.0) * dCdt2, dCdt3);
	derivative(C_CIS +  dt        * dCdt3, dCdt4);

	C_CIS += (dt / 6.0) * (dCdt1 + 2.0 * dCdt2 + 2.0 * dCdt3 + dCdt4);
    }

    void run()
    {
	cout << "Start quantum simulation" << endl;

	size_t ncount = 0, nstepmax = size_t(tmax/dt+1);
	for(double t = 0; t < tmax; t+=dt)
	{
	    Observe(t, ncount, nstepmax);

	    // propagate in the time domain
	    //Step_RK4();
	    C_CIS = expHCISdt * C_CIS;
	    ncount ++;
	}

	// finally output results
	traj_popu.save("traj_QM.txt", arma::raw_ascii);
	traj_Ez.save("Ez_QM.txt", arma::raw_ascii);
	traj_Ez2.save("Ez2_QM.txt", arma::raw_ascii);
	
	cout << "END quantum simulation" << endl;
    }
};



int main()
{
    QuantumSimulator model_QM;
    model_QM.run();

    return 0;
}
