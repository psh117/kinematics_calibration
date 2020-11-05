
#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <mutex>
#include <thread>
#include <robot_model/franka_panda_model.h>
#include "suhan_benchmark.h"

class HCTest
{
  const int N_DH = 4; // number of dh parameters to calibrate
  const int N_J = 7; // number of joints to calibrate
  const int N_CAL = N_DH*N_J; // total variables to calibrate
  const double d_lr = 0.15;
  const double d_lt = 0.15;
  const double d_rt = 0.15;

  struct state_
  {
    Eigen::VectorXd alpha(N_J)
    Eigen::VectorXd a(N_J)
    Eigen::VectorXd d(N_J);
    Eigen::VectorXd offset(N_CAL); // x(a,d,theta,alpha)
    Eigen::VectorXd p_total;
    Eigen::MatrixXd jac_total;
    Eigen::VectorXd del_phi;
  };

  std::map<std::string, std::shared_ptr<a_state_>> arm_state;
  state_ left_state_;
  state_ right_state_;
  state_ top_state_;

  std::mutex data_mutex;
  std::mutex print_mutex;


public:
  Eigen::Isometry3d transformFromDH(double a, double d, double alpha, double theta);

private:
  
};