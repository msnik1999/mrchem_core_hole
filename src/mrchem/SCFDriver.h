#ifndef SCFDRIVER_H
#define SCFDRIVER_H

#include <vector>
#include <string>
#include <Eigen/Core>

class Getkw;

class Nuclei;
class Molecule;
class OrbitalSet;
class GroundStateSolver;
class LinearResponseSolver;
class HelmholtzOperatorSet;
class Accelerator;

class QMOperator;
class FockOperator;
class PoissonOperator;
class CoulombPotential;
class CoulombHessian;
class KineticOperator;
class NuclearPotential;
class ExchangePotential;
class ExchangeHessian;
class XCPotential;
class XCHessian;
class XCFunctional;

class SCFDriver {
public:
    SCFDriver(Getkw &input);
    virtual ~SCFDriver() { }

    void setup();
    void run();
    void clear();

protected:
    double est_norm;
    double rel_prec;
    double r_O[3];
    bool center_of_mass;
    std::vector<double> gauge;

    // Run parameters
    bool run_el_field_rsp;
    bool run_mag_field_rsp;
    bool run_mag_moment_rsp;
    bool run_ground_state;
    bool run_dipole_moment;
    bool run_quadrupole_moment;
    bool run_polarizability;
    bool run_optrot_electric;
    bool run_optrot_magnetic;
    bool run_magnetizability;
    bool run_nmr_shielding;
    bool run_spin_spin;
    bool velocity_gauge;
    std::vector<int> nmr_nuclei;
    std::vector<int> spin_spin_k;
    std::vector<int> spin_spin_l;
    std::vector<double> frequencies;

    // Molecule input
    int mol_charge;
    int mol_multiplicity;
    std::vector<std::string> mol_coords;

    // Wavefunction input
    bool wf_restricted;
    std::string wf_method;

    // DFT input
    bool dft_spin;
    double dft_x_fac;
    std::vector<double> dft_cutoff;
    std::vector<double> dft_func_coefs;
    std::vector<std::string> dft_func_names;

    // Ground state input
    std::string scf_start;
    std::string scf_acc;
    int scf_history;
    int scf_max_iter;
    int scf_rotation;
    bool scf_localize;
    bool scf_write_orbitals;
    double scf_orbital_thrs;
    double scf_property_thrs;
    double scf_lambda_thrs;
    std::vector<double> scf_orbital_prec;

    // Response input
    std::string rsp_start;
    std::string rsp_acc;
    int rsp_history;
    int rsp_max_iter;
    bool rsp_localize;
    bool rsp_write_orbitals;
    double rsp_orbital_thrs;
    double rsp_property_thrs;
    std::vector<int> rsp_dir;
    std::vector<double> rsp_orbital_prec;

    // File input
    std::string file_start_orbitals;
    std::string file_final_orbitals;
    std::string file_start_x_orbs;
    std::string file_final_x_orbs;
    std::string file_start_y_orbs;
    std::string file_final_y_orbs;
    std::string file_basis_set;
    std::string file_dens_mat;
    std::string file_fock_mat;
    std::string file_energy_vec;
    std::string file_mo_mat_a;
    std::string file_mo_mat_b;

    // SCF machinery
    HelmholtzOperatorSet *helmholtz;
    Accelerator *scf_kain;
    Accelerator *rsp_kain_x;
    Accelerator *rsp_kain_y;

    // Unperturbed quantities
    Molecule *molecule;
    Nuclei *nuclei;
    OrbitalSet *orbitals;
    PoissonOperator *P;
    KineticOperator *T;
    NuclearPotential *V;
    CoulombPotential *J;
    ExchangePotential *K;
    XCPotential *XC;
    XCFunctional *xcfun_1;
    FockOperator *f_oper;
    Eigen::MatrixXd *f_mat;

    // Perturbed quantities
    OrbitalSet *x_orbs;
    OrbitalSet *y_orbs;
    CoulombHessian *dJ;
    ExchangeHessian *dK;
    XCHessian *dXC;
    XCFunctional *xcfun_2;
    FockOperator *df_oper;

    bool sanityCheck() const;
    bool runInitialGuess(FockOperator &oper, Eigen::MatrixXd &F, OrbitalSet &orbs);
    bool runGroundState();
    void runElectricFieldResponse(double omega);
    void runMagneticFieldResponse(double omega);
    void runMagneticMomentResponse(const std::string &type, int L);

    void calcGroundStateProperties();
    void calcElectricFieldProperties(int d, double omega);
    void calcMagneticFieldProperties(int d, double omega);
    void calcMagneticMomentProperties(const std::string &type, int d, int L);

    GroundStateSolver *setupInitialGuessSolver();
    GroundStateSolver *setupGroundStateSolver();
    LinearResponseSolver *setupLinearResponseSolver(bool dynamic);

    void setupInitialGroundState();
    void setupInitialResponse(QMOperator &h, int d, bool dynamic, double omega);

    void setupPerturbedOrbitals(bool dynamic);
    void clearPerturbedOrbitals(bool dynamic);

    void setupPerturbedOperators(QMOperator &dH, bool dynamic);
    void clearPerturbedOperators();
};

#endif
