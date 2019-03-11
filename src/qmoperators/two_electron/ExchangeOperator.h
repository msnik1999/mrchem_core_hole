#pragma once

#include "ExchangePotential.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class ExchangeOperator
 *
 * @brief Operator containing a single ExchangePotential
 *
 * The operator is defined as the Hartree-Fock exchange arising from a particular
 * set of orbitals. The OrbitalVector defining the operator is fixed throughout the
 * operators life time, but the orbitals themselves are allowed to change in between
 * each application. The internal exchange potentials (the operator applied to it's
 * own orbitals) can be precomputed and stored for fast retrieval.
 */

namespace mrchem {

class ExchangeOperator final : public RankZeroTensorOperator {
public:
    ExchangeOperator(std::shared_ptr<mrcpp::PoissonOperator> P, std::shared_ptr<OrbitalVector> Phi)
            : exchange(new ExchangePotential(P, Phi, false)) {
        RankZeroTensorOperator &K = (*this);
        K = exchange;
    }
    ~ExchangeOperator() override = default;

    auto &getPoisson() { return exchange->getPoisson(); }
    void setupInternal(double prec) { exchange->setupInternal(prec); }
    void rotate(const ComplexMatrix &U) { exchange->rotate(U); }

    ComplexDouble trace(OrbitalVector &Phi) { return 0.5 * RankZeroTensorOperator::trace(Phi); }

protected:
    std::shared_ptr<ExchangePotential> exchange;
};

} // namespace mrchem
