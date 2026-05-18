/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "qmoperator_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/QMPotential.h"
#include "qmoperators/one_electron/MomentumOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace qmoperator {
//Non-relativistic
ComplexMatrix calc_kinetic_matrix_component(int d, MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket);
//scalar relativistic 
ComplexMatrix calc_kinetic_matrix_component(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket);
//spinorial relativistic
ComplexMatrix calc_kinetic_matrix_component(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket);
//scalar relativistic
ComplexMatrix calc_kinetic_matrix_component_symmetrized(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial = false);
} // namespace qmoperator

double qmoperator::calc_kinetic_trace(MomentumOperator &p, OrbitalVector &Phi) {
    DoubleVector eta = orbital::get_occupations(Phi).cast<double>();
    DoubleVector norms = DoubleVector::Zero(Phi.size());
    {
        OrbitalVector dPhi = p[0](Phi);
        norms += orbital::get_squared_norms(dPhi);
    }
    {
        OrbitalVector dPhi = p[1](Phi);
        norms += orbital::get_squared_norms(dPhi);
    }
    {
        OrbitalVector dPhi = p[2](Phi);
        norms += orbital::get_squared_norms(dPhi);
    }
    return 0.5 * eta.dot(norms);
}

ComplexDouble qmoperator::calc_kinetic_trace(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &Phi, bool spinorial) {
    ComplexDouble out = {0.0, 0.0};
    // MSG_INFO("kin trace start");
    int alpha_index = 0; // for spinorial operators, the kinetic operator is of the form sigma p V sigma p, with sigma being a Pauli matrix. alpha_index represents which Pauli matrix is used (x,y or z). For non-spinorial operators, it is unused.
    {   
        if (spinorial) alpha_index = 1; // if spinorial, we use sigma_x for the kinetic operator, but it is arbitrary which Pauli matrix we use, as long as it is the same on the left and right of V.
        OrbitalVector dPhi = p[0](Phi, alpha_index);
        out += V.trace(dPhi);
    }
    {   
        if (spinorial) alpha_index = 2;
        OrbitalVector dPhi = p[1](Phi, alpha_index);
        out += V.trace(dPhi);
    }
    {
        if (spinorial) alpha_index = 3;
        OrbitalVector dPhi = p[2](Phi, alpha_index);
        out += V.trace(dPhi);
    }
    return 0.5 * out;
}

/** @brief Expectation value matrix: T_ij = <i|T|j> = <i|p p|j>
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix qmoperator::calc_kinetic_matrix(MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket) { 
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component(0, p, bra, ket);
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component(1, p, bra, ket);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component(2, p, bra, ket);
    return T_x + T_y + T_z;
}

/** @brief Expectation value matrix ZORA: T_ij = <i|T_zora|j> = <i|p V_zora p|j>
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix qmoperator::calc_kinetic_matrix(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    if (spinorial and (bra[0].Ncomp() < 2 or ket[0].Ncomp() < 2)) MSG_WARN("Spin-orbit interaction enabled with scalar functions, will probably create issues");
    MSG_INFO("Spinorial = "<< spinorial);
    if (not spinorial){
        //scalar case
        ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component(0, p, V, bra, ket); 
        ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component(1, p, V, bra, ket);
        ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component(2, p, V, bra, ket);
        return T_x + T_y + T_z;
    } else {
        //spinorial case
        return qmoperator::calc_kinetic_matrix_component(p, V, bra, ket);
    }
}

ComplexMatrix qmoperator::calc_kinetic_matrix_symmetrized(MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    // MSG_INFO("calc kin mat symm");
    ComplexMatrix T_x = qmoperator::calc_kinetic_matrix_component_symmetrized(0, p, V, bra, ket, spinorial);
    ComplexMatrix T_y = qmoperator::calc_kinetic_matrix_component_symmetrized(1, p, V, bra, ket, spinorial);
    ComplexMatrix T_z = qmoperator::calc_kinetic_matrix_component_symmetrized(2, p, V, bra, ket, spinorial);
    return T_x + T_y + T_z;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component(int d, MomentumOperator &p, OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);
    // MSG_INFO("calc kin mat comp start braket true=" << (&bra == &ket));

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = p[d](ket);
        // MSG_INFO("1");
        nNodes += orbital::get_n_nodes(dKet);
        // MSG_INFO("2");
        sNodes += orbital::get_size_nodes(dKet);
        // MSG_INFO("3");
        T = mrcpp::calc_overlap_matrix(dKet);
    } else {
        OrbitalVector dBra = p[d](bra);
        OrbitalVector dKet = p[d](ket);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = mrcpp::calc_overlap_matrix(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

ComplexMatrix qmoperator::calc_kinetic_matrix_component_symmetrized(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket, bool spinorial) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    //in case of spinorial orbitals, the kinetic operator is sigma dot p, we therefore need to set the index of the pauli matrix to use
    int pauli_index = 0; //identity (fits non-spinorial case)
    if (spinorial) pauli_index = d + 1; //Pauli index is 1 for x, 2 for y and 3 for z, whereas the momentum operator index is 0 for x, 1 for y and 2 for z, so we need to add 1 to get the correct Pauli matrix.

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        OrbitalVector dKet = (V * p[d])(ket, pauli_index);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = mrcpp::calc_overlap_matrix(dKet, dKet);
    } else {
        OrbitalVector dBra = (V * p[d])(bra, pauli_index);
        OrbitalVector dKet = (V * p[d])(ket, pauli_index);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = mrcpp::calc_overlap_matrix(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

/* @brief Scalar relativistic computation of the kinetic matrix
* @param d: integer for selecting the direction, x,y,z
*/
ComplexMatrix qmoperator::calc_kinetic_matrix_component(int d, MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);

    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        // scalar case, so no need to apply pauli matrices (second argument, index=0 means applying the identity matrix)
        OrbitalVector dKet = p[d](ket, 0); 

        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        T = V(dKet, dKet); 
    } else {
        // scalar case, so no need to apply pauli matrices (second argument, index=0 means applying the identity matrix)
        OrbitalVector dBra = p[d](bra, 0);
        OrbitalVector dKet = p[d](ket, 0);
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        T = V(dBra, dKet);
    }
    if (d == 0) mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 1) mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    if (d == 2) mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    return 0.5 * T;
}

/* @brief Spinorial computation of the kinetic matrix 
* The directions can no longer be computed seperately, due to the kinetic operator becoming σ·p
*/
ComplexMatrix qmoperator::calc_kinetic_matrix_component( MomentumOperator &p, RankZeroOperator &V, OrbitalVector &bra, OrbitalVector &ket) {
    Timer timer;
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T = ComplexMatrix::Zero(Ni, Nj);
    // MSG_INFO("calc kin mat comp start d=" << d << " Ni=" << Ni << " Nj=" << Nj);

    // MSG_INFO("spinorial done");
    int nNodes = 0, sNodes = 0;
    if (&bra == &ket) {
        // MSG_INFO("bra ket same");
        OrbitalVector dKet_x = p[0](ket, 1); //(σ_x·p_x)|ket>
        OrbitalVector dKet_y = p[1](ket, 2); //(σ_y·p_y)|ket>
        OrbitalVector dKet_z = p[2](ket, 3); //(σ_z·p_z)|ket>
        //summing it all together
        OrbitalVector dKet_tmp = orbital::add({1.0,0.0}, dKet_x, {1.0,0.0}, dKet_y);//intermediate sum
        OrbitalVector dKet = orbital::add({1.0,0.0}, dKet_tmp, {1.0,0.0}, dKet_z);
        // ComplexVector ones = ComplexVector::Ones(Ni, 1); //Vector of 1.0 for in-place addition later
        // OrbitalVector dKet = p[1](ket, 1); //(σ_x·p_x)|ket>, also holding the sum of all three directions
        // // mrcpp::add(dKet, ones, p[2](ket, 2)); //(σ_y·p_y)|ket>
        // // mrcpp::add(dKet, ones, p[3](ket, 3)); //(σ_z·p_z)|ket>
        // OrbitalVector dKet_y = p[2](ket, 2); //(σ_y·p_y)|ket>
        // mrcpp::add(dKet, ones, dKet_y); 
        // OrbitalVector dKet_z = p[3](ket, 3); //(σ_z·p_z)|ket>
        // mrcpp::add(dKet, ones, dKet_z);
        
        // ComplexMatrix testKet = mrcpp::calc_overlap_matrix(dKet);
        // // MSG_INFO("bk same dket test " << testKet(0,0) << " "<< testKet(0,1) << " " << d);
        // // MSG_INFO("bk same dket test " << testKet(1,0) << " "<< testKet(1,1) << " " << d);
        // OrbitalVector dKettest = p[d](ket, 0); //debug test
        // ComplexMatrix testKet2 = mrcpp::calc_overlap_matrix(dKet); //debug test
        // // MSG_INFO("bk same dket test 2" << testKet2(0,0) << " "<< testKet2(0,1) << " " << d);
        // // MSG_INFO("bk same dket test 2" << testKet2(1,0) << " "<< testKet2(1,1) << " " << d);
        // //further debug stuff
        // ComplexMatrix Ttest = ComplexMatrix::Zero(Ni,Nj);
        // Ttest = V(dKettest, dKettest);
        // // MSG_INFO("bk same V applied test "  << Ttest(0,0) << " "<< Ttest(0,1) << " " << d);
        // // MSG_INFO("bk same V applied test "  << Ttest(1,0) << " "<< Ttest(1,1));
        // //end further debug sutff

        // MSG_INFO("bk same dket done ");
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dKet);
        // MSG_INFO("bk same gotten stuff ");
        T = V(dKet, dKet); //includes conjugation?
        // MSG_INFO("test")
        // MSG_INFO("bk same V applied tut "  << T(0,0) << " "<< T(0,1) << " " << d);
        // MSG_INFO("bk same V applied tut "  << T(1,0) << " "<< T(1,1));
        
    } else {
        // MSG_INFO("bra no ket");
        OrbitalVector dBra_x = p[0](bra, 1); //(σ_x·p_x)|ket>
        OrbitalVector dBra_y = p[1](bra, 2); //(σ_y·p_y)|ket>
        OrbitalVector dBra_z = p[2](bra, 3); //(σ_z·p_z)|ket>
        //summing it all together
        OrbitalVector dBra_tmp = orbital::add({1.0,0.0}, dBra_x, {1.0,0.0}, dBra_y);//intermediate sum
        OrbitalVector dBra = orbital::add({1.0,0.0}, dBra_tmp, {1.0,0.0}, dBra_z);

        OrbitalVector dKet_x = p[0](ket, 1); //<bra|(σ_x·p_x)
        OrbitalVector dKet_y = p[1](ket, 2); //<bra|(σ_y·p_y)
        OrbitalVector dKet_z = p[2](ket, 3); //<bra|(σ_z·p_z)
        //summing it all together
        OrbitalVector dKet_tmp = orbital::add({1.0,0.0}, dKet_x, {1.0,0.0}, dKet_y);//intermediate sum
        OrbitalVector dKet = orbital::add({1.0,0.0}, dKet_tmp, {1.0,0.0}, dKet_z);
        // ComplexVector ones = ComplexVector::Ones(Ni, 1);
        // OrbitalVector dBra = p[1](bra, 1); //(σ_x·p_x)|ket>, also holding the sum of all three directions
        // // mrcpp::add(dBra, ones, p[2](bra, 2)); //(σ_y·p_y)|ket>
        // // mrcpp::add(dBra, ones, p[3](bra, 3)); //(σ_z·p_z)|ket>
        // OrbitalVector dBra_y = p[2](bra, 2); //(σ_y·p_y)|ket>
        // mrcpp::add(dBra, ones, dBra_y); 
        // OrbitalVector dBra_z = p[3](bra, 3); //(σ_z·p_z)|ket>
        // mrcpp::add(dBra, ones, dBra_z);
        
        // ones = ComplexVector::Ones(Nj, 1);
        // OrbitalVector dKet = p[1](ket, 1); //(σ_x·p_x)|ket>, also holding the sum of all three directions
        // // mrcpp::add(dKet, ones, p[2](ket, 2)); //(σ_y·p_y)|ket>
        // // mrcpp::add(dKet, ones, p[3](ket, 3)); //(σ_z·p_z)|ket>
        // OrbitalVector dKet_y = p[2](ket, 2); //(σ_y·p_y)|ket>
        // mrcpp::add(dKet, ones, dKet_y); 
        // OrbitalVector dKet_z = p[3](ket, 3); //(σ_z·p_z)|ket>
        // mrcpp::add(dKet, ones, dKet_z);
        // MSG_INFO("bk diff ket done");
        nNodes += orbital::get_n_nodes(dBra);
        nNodes += orbital::get_n_nodes(dKet);
        sNodes += orbital::get_size_nodes(dBra);
        sNodes += orbital::get_size_nodes(dKet);
        // MSG_INFO("bk diff gotten stuff");
        T = V(dBra, dKet);
        // MSG_INFO("bk diff V applied ");
    }
    // MSG_INFO("calc kin mat out of if statement");
    mrcpp::print::tree(2, "<i|sigma p kappa sigma p|j>", nNodes, sNodes, timer.elapsed());
    // MSG_INFO("end ");
    return 0.5 * T;
}

} // namespace mrchem
