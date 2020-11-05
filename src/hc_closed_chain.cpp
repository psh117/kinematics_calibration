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

// #define TEST_PRINT

const int N_PARAM = 3;
const int N_ARM = 3;  //number of arms
const int N_DH = 4;   // number of dh parameters to calibrate
const int N_J = 7;    // number of joints in one robot
const int N_CAL = N_J * N_DH * N_ARM; // total variables to calibrate
const int N_JDH = N_J * N_DH;    // number of joints in all robots

// each vector's information is each arm's information
std::vector<Eigen::Matrix<double, N_J, N_DH>> offset_matrix;  // each col (a, d, theta, alpha)
std::vector<Eigen::Matrix<double, N_ARM, N_J>> theta_data;    // dataset
std::vector<Eigen::Isometry3d> T_W0;
std::vector<FrankaPandaModel> fpm;
std::vector<double> distTrue;

std::string data_name = "data_1";
Eigen::VectorXd del_phi(N_CAL);
Eigen::VectorXd p_total;
Eigen::MatrixXd jacobian;
int num_data;

void initialize()
{
  std::ifstream rf; 
  rf.open(data_name+".txt");
  while (!rf.eof())
  {  
    Eigen::Matrix<double, N_ARM, N_J> ar;
    for (int i=0; i<N_ARM; i++)
    {
      Eigen::Matrix<double, 1, N_J> d;
      for (int j=0; j<N_J; j++)
      {
        rf >> d(j);
      }
      ar.row(i) = d;
    }
    theta_data.push_back(ar);
  }
  rf.close();

  num_data = theta_data.size();
  jacobian.resize(N_PARAM*num_data, N_CAL);
  p_total.resize(N_PARAM*num_data);

  Eigen::Isometry3d T_W0_LEFT;
  Eigen::Isometry3d T_W0_RIGHT;
  Eigen::Isometry3d T_W0_TOP;
  T_W0_LEFT.linear() << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;
  T_W0_LEFT.translation() << 0.0, 0.3, 1.0;
  T_W0_RIGHT.linear() << 1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0,
                        0.0, 0.0, 1.0;
  T_W0_RIGHT.translation() << 0.0, -0.3, 1.0;
  T_W0_TOP.linear() << -1.0,  0.0, 0.0,
                        0.0, -1.0, 0.0,
                        0.0,  0.0, 1.0;
  T_W0_TOP.translation() << 1.35, 0.3, 1.0;
  T_W0.push_back(T_W0_LEFT);
  T_W0.push_back(T_W0_RIGHT);
  T_W0.push_back(T_W0_TOP);

  distTrue.push_back(0.6); // LEFT&RIGHT
  distTrue.push_back(0.736078); // RIGHT&TOP
  distTrue.push_back(0.426392); // TOP&LEFT

  for (int i=0; i<N_ARM; i++)
  {
    FrankaPandaModel fpm_ = FrankaPandaModel();
    Eigen::Matrix<double, N_J, N_DH> x;
    x.setZero();

    offset_matrix.push_back(x);
    fpm.push_back(fpm_);
  }

  std::cout << "Number of datasets: " << theta_data.size() << std::endl;
  std::cout << "Distance between L&R: " << distTrue[0] << std::endl;
  std::cout << "Distance between R&T: " << distTrue[1] << std::endl;
  std::cout << "Distance between T&L: " << distTrue[2] << std::endl;
}

Eigen::VectorXd getDistanceDiff()
{
  Eigen::VectorXd del_(N_PARAM*num_data);
  del_.setZero();

  for (int i=0; i<num_data; i++)
  {
    auto T0 = (T_W0[0]) * fpm[0].getTransform(theta_data[i].row(0));
    auto T1 = (T_W0[1]) * fpm[0].getTransform(theta_data[i].row(1));
    auto T2 = (T_W0[2]) * fpm[0].getTransform(theta_data[i].row(2));

    Eigen::Vector3d X0 = T0.translation();
    Eigen::Vector3d X1 = T1.translation();
    Eigen::Vector3d X2 = T2.translation();

    del_(i*N_PARAM + 0) = distTrue[0] - (X0 - X1).norm();
    del_(i*N_PARAM + 1) = distTrue[1] - (X1 - X2).norm();
    del_(i*N_PARAM + 2) = distTrue[2] - (X2 - X0).norm();

    if (N_PARAM > 3)
    {
      Eigen::Quaterniond q0(T0.linear());
      Eigen::Quaterniond q1(T1.linear());
      Eigen::Quaterniond q2(T2.linear() * T_W0[2].linear());

      del_(i*N_PARAM + 3) = q0.angularDistance(q1);
      del_(i*N_PARAM + 4) = q1.angularDistance(q2);
      del_(i*N_PARAM + 5) = q2.angularDistance(q0);
    }

#ifdef TEST_PRINT
    std::cout<<"Transform_LEFT: "<<X0.transpose()<<std::endl;
    std::cout<<"Transform_RIGHT: "<<X1.transpose()<<std::endl;
    std::cout<<"Transform_TOP: "<<X2.transpose()<<std::endl;
    std::cout<<"(X0 - X1).norm(): "<<(X0 - X1).norm()<<std::endl;
    std::cout<<"(X1 - X2).norm(): "<<(X1 - X2).norm()<<std::endl;
    std::cout<<"(X2 - X0).norm(): "<<(X2 - X0).norm()<<std::endl;
    std::cout<<"del_(i*N_PARAM + 0): "<<del_(i*N_PARAM + 0)<<std::endl;
    std::cout<<"del_(i*N_PARAM + 1): "<<del_(i*N_PARAM + 1)<<std::endl;
    std::cout<<"del_(i*N_PARAM + 2): "<<del_(i*N_PARAM + 2)<<std::endl;
    if (N_PARAM > 3)
    {
      std::cout<<"del_(i*N_PARAM + 3): "<<del_(i*N_PARAM + 3)<<std::endl;
      std::cout<<"del_(i*N_PARAM + 4): "<<del_(i*N_PARAM + 4)<<std::endl;
      std::cout<<"del_(i*N_PARAM + 5): "<<del_(i*N_PARAM + 5)<<std::endl;
    }
#endif
  }
  return del_;
}

void getJacobian()
{
  Eigen::VectorXd t1, t2, m1, m2, m3;
  for (int arm_=0; arm_<N_ARM; arm_++)
  {
    Eigen::Matrix<double, N_J, N_DH> y1 = offset_matrix[arm_];
    Eigen::Matrix<double, N_J, N_DH> y2 = offset_matrix[arm_];
    for (int j=0; j<N_J; j++)
    {
      for (int dh_=0; dh_<N_DH; dh_++)
      {
        const double ax = std::fabs(offset_matrix[arm_](j, dh_));
        // Make step size as small as possible while still giving usable accuracy.
        const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);
        y1(j, dh_) += h;
        fpm[arm_].initModel(y1);
        t1 = getDistanceDiff();
        y2(j, dh_) -= h;
        fpm[arm_].initModel(y2);
        t2 = getDistanceDiff();
        m1 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        y1(j, dh_) += h;
        fpm[arm_].initModel(y1);
        t1 = getDistanceDiff();
        y2(j, dh_) -= h;
        fpm[arm_].initModel(y2);
        t2 = getDistanceDiff();
        m2 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        y1(j, dh_) += h;
        fpm[arm_].initModel(y1);
        t1 = getDistanceDiff();
        y2(j, dh_) -= h;
        fpm[arm_].initModel(y2);
        t2 = getDistanceDiff();
        m3 = (t1 - t2) / (y1(j, dh_) - y2(j, dh_));

        jacobian.col(arm_*N_JDH + j*N_DH + dh_) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

        // Reset for next iteration.
        y1(j, dh_) = y2(j, dh_) = offset_matrix[arm_](j, dh_);
      }
    }
    // Reset for next iteration.
    fpm[arm_].initModel(offset_matrix[arm_]);
  }
}

int main(int argc, char**argv)
{
  initialize();

  int iter = 100;
  while (iter--)
  {
    for (int i=0; i<N_ARM; i++) {fpm[i].initModel(offset_matrix[i]);}
    p_total = getDistanceDiff();
    getJacobian();
    std::cout << "\neval: " << p_total.squaredNorm() / num_data << std::endl;
    del_phi = jacobian.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p_total);
    std::cout << "\ndel_phi.norm(): " << del_phi.norm() << std::endl;
    for (int arm_=0; arm_<N_ARM; arm_++)
    {
      for (int j=0; j<N_J; j++)
      {
        offset_matrix[arm_].row(j) -= del_phi.segment<N_DH>(arm_*N_JDH + j*N_DH);
      }
    }
    if (del_phi.norm() < 1e-9) 
    {
      std::cout << "\noffset_matrix LEFT:\n" << offset_matrix[0] << std::endl;
      std::cout << "\noffset_matrix RIGHT:\n" << offset_matrix[1] << std::endl;
      std::cout << "\noffset_matrix TOP:\n" << offset_matrix[2] << std::endl;
      std::cout<<"reached optimal value at iter: "<<100 - iter<<std::endl;
      break;
    }
  }
  return 0;
}