#ifndef GROUNDSTATESOLVER_H
#define GROUNDSTATESOLVER_H

#include <Eigen/Core>

#include "SCF.h"
#include "SCFEnergy.h"

class Accelerator;

class GroundStateSolver : public SCF {
public:
    GroundStateSolver(HelmholtzOperatorSet &h, Accelerator *a = 0);
    virtual ~GroundStateSolver();

    void setup(FockOperator &f_oper, Eigen::MatrixXd &f_mat, OrbitalVector &phi);
    void clear();

    bool optimize();

protected:
    std::vector<SCFEnergy> energy;

    FockOperator *fOper_n;
    FockOperator *fOper_np1;

    OrbitalVector *orbitals_n;
    OrbitalVector *orbitals_np1;
    OrbitalVector *dOrbitals_n;

    Eigen::MatrixXd *fMat_n;
    Eigen::MatrixXd *fMat_np1;
    Eigen::MatrixXd *dfMat_n;

    Accelerator *accelerator;	    // Pointer to external object, do not delete!

    bool optimizeOrbitals();
    bool optimizeEnergy();

    void setupUpdates();
    void clearUpdates();

    void rotate(OrbitalVector &phi,
                Eigen::MatrixXd &f_mat,
                FockOperator *f_oper = 0);

    void calcOrbitalUpdates();
    void calcFockMatrixUpdate();

    Orbital* getHelmholtzArgument(int i,
                                  OrbitalVector &phi,
                                  Eigen::MatrixXd &f_mat,
                                  bool adjoint);

    void printProperty() const;
    int printTreeSizes() const;

    double calcProperty();
    double calcTotalError() const;
    double calcOrbitalError() const;
    double calcPropertyError() const;
};

#endif // GROUNDSTATESOLVER_H
