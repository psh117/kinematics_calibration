
#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <mutex>
#include <thread>
#include <robot_model/franka_panda_model.h>

const int N_PARAM = 6;
const int N_ARM = 3;  //number of arms
const int N_DH = 4;   // number of dh parameters to calibrate
const int N_J = 7;    // number of joints in one robot
const int N_CAL = N_J * N_DH * N_ARM; // total variables to calibrate
const int N_JDH = N_J * N_DH;    // number of joints in all robots

struct state_
{
  Eigen::VectorXd pose_total;
  Eigen::VectorXd del_pose_total;
  Eigen::VectorXd del_phi;
  Eigen::Matrix<double, N_J, 4> offset_matrix; // each col (a, d, theta, alpha)
  Eigen::MatrixXd jacobian;
  Eigen::Isometry3d T_W0;
  std::vector<Eigen::Matrix<double, N_J, 1>> theta;
  std::string arm_name;
  std::string arm_true;
  FrankaPandaModel fpm;
};

Eigen::Matrix<double, N_J, N_DH> original_dh_matrix;  // each col (a, d, theta, alpha)
Eigen::Matrix<double, N_JDH, N_JDH> weight;
int num_data;
double EE = 0.107;

//   /home/kimhc/git/kinematics_calibration/data/CLOSED_CALIBRATION/
std::string current_workspace = "/home/kimhc/git/kinematics_calibration/data/triple_calibration/";
std::string data_input;
std::string data_iter;
std::string data_offset;
std::string offset_loader;
Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");

std::map<std::string, std::shared_ptr<state_>> arm_state;
state_ left_state_;
state_ right_state_;
state_ top_state_;