#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NuclearOperator.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief projects an analytic expression for a smoothed nuclear potential
 *
 * @param[in] nucs collection of nuclei that defines the potential
 * @param[in] prec precision used both in smoothing and projection
 *
 * We create two different analytic functions for the nuclear potential:
 *
 * 1) this->func: The total potential from all nuclei of the system.
 *                This is needed later for analytic calculations.
 *
 * 2) loc_func: Temporary function that contains only some of the
 *              nuclei, which are distributed among the available
 *              MPIs. This is used only for the projection below.
 */
NuclearPotential::NuclearPotential(const Nuclei &nucs, double proj_prec, double smooth_prec, bool mpi_share)
        : QMPotential(1, mpi_share) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;

    int oldprec = Printer::setPrecision(5);
    mrcpp::print::header(0, "Building nuclear potential");
    println(0, "    N   Atom            Charge      Precision    Smoothing ");
    mrcpp::print::separator(0, '-');

    NuclearFunction loc_func;

    double c = 0.00435 * smooth_prec;
    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double Z_5 = std::pow(Z, 5.0);
        double smooth = std::pow(c / Z_5, 1.0 / 3.0);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;

        this->func.push_back(nuc, smooth);
        if (mpi::orb_rank == proj_rank) loc_func.push_back(nuc, smooth);

        std::stringstream o_sym;
        o_sym << std::setw(4) << k;
        o_sym << std::setw(6) << nuc.getElement().getSymbol();
        std::array<double, 3> array{Z, smooth_prec, smooth};
        print_utils::coord(0, o_sym.str(), array, 5, true);
    }

    Timer t_tot;

    // Scale precision by system size
    int Z_tot = chemistry::get_total_charge(nucs);
    double abs_prec = proj_prec / Z_tot;

    QMFunction V_loc(false);

    Timer t_loc;
    qmfunction::project(V_loc, loc_func, NUMBER::Real, abs_prec);
    t_loc.stop();

    Timer t_com;
    allreducePotential(abs_prec, V_loc);
    t_com.stop();

    t_tot.stop();
    mrcpp::print::separator(0, '-');
    print_utils::qmfunction(0, "Local potential", V_loc, t_loc);
    print_utils::qmfunction(0, "Allreduce", V_loc, t_com);
    mrcpp::print::footer(0, t_tot, 2);
    Printer::setPrecision(oldprec);
}

void NuclearPotential::allreducePotential(double prec, QMFunction &V_loc) {
    QMFunction &V_tot = *this;

    // Add up local contributions into the grand master
    mpi::reduce_function(prec, V_loc, mpi::comm_orb);
    if (mpi::grand_master()) {
        // If numerically exact the grid is huge at this point
        if (mpi::numerically_exact) V_loc.crop(prec);
    }

    if (not V_tot.hasReal()) V_tot.alloc(NUMBER::Real);
    if (V_tot.isShared()) {
        int tag = 3141;
        // MPI grand master distributes to shared masters
        mpi::broadcast_function(V_loc, mpi::comm_sh_group);
        if (mpi::share_master()) {
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(V_tot.real(), V_loc.real());
            mrcpp::copy_func(V_tot.real(), V_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mpi::share_function(V_tot, 0, tag, mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mpi::broadcast_function(V_loc, mpi::comm_orb);
        // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(V_tot.real(), V_loc.real());
        mrcpp::copy_func(V_tot.real(), V_loc.real());
    }
}

/** @brief computes the interaction energy of the nuclear potential with a second set of nuclei.
 *
 * @param[in] nucs the set of nuclei
 *
 * Note: this function is not suited to compute the nuclear self-energy
 *
 */
double NuclearOperator::trace(const Nuclei &nucs) {
    MSG_WARN("This routine has never been tested!");
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int k = 0; k < nNucs; k++) {
        const Nucleus &nuc_k = nucs[k];
        double Z_k = nuc_k.getCharge();
        const mrcpp::Coord<3> &R_k = nuc_k.getCoord();
        E_nuc += Z_k * this->r_m1->evalf(R_k);
    }
    return E_nuc;
}

} // namespace mrchem
