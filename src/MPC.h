#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  // This is the length from front to CoG that has a similar radius.
  double Lf_;
  // TODO (DONE): Set the timestep length and duration
  size_t N_;
  double dt_;
  
  vector<double> Xpred_;
  vector<double> Ypred_;
  
  MPC();

  virtual ~MPC();
  
  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
