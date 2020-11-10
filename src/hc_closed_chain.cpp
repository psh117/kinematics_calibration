#include <iostream>
#include <fstream>
#include <chrono>
#include <Eigen/Dense>
#include <map>
#include <vector>
#include <mutex>
#include <thread>
#include <string>
#include <robot_model/franka_panda_model.h>
#include "suhan_benchmark.h"

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

// // Generic functor
// template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
// struct Functor
// {
// typedef _Scalar Scalar;
// enum {
//     InputsAtCompileTime = NX,
//     ValuesAtCompileTime = NY
// };
// typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
// typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
// typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

// int m_inputs, m_values;

// Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
// Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

// int inputs() const { return m_inputs; }
// int values() const { return m_values; }

// };

// struct my_functor : Functor<double>
// {
// my_functor(void): Functor<double>(2,2) {}
// int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
// {
//     // Implement y = 10*(x0+3)^2 + (x1-5)^2
//     fvec(0) = 10.0*pow(x(0)+3.0,2) +  pow(x(1)-5.0,2);
//     fvec(1) = 0;

//     return 0;
// }
// };


// #define DEBUG_MODE

const int N_PARAM = 6;
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

//   /home/kimhc/git/kinematics_calibration/data/CLOSED_CALIBRATION/
std::string user_name = "kimhc";
std::string current_workspace = "/home/" + user_name + "/git/kinematics_calibration/";
std::string data_input;
std::string data_iter;
std::string data_offset;
Eigen::VectorXd del_phi(N_CAL);
Eigen::VectorXd p_total;
Eigen::MatrixXd jacobian;
int num_data;
Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");

void initialize()
{
  std::ifstream rf; 
  rf.open(data_input);
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

  distTrue.push_back(0.207); // LEFT&RIGHT
  distTrue.push_back(0.292742); // RIGHT&TOP
  distTrue.push_back(0.207); // TOP&LEFT

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
    auto T1 = (T_W0[1]) * fpm[1].getTransform(theta_data[i].row(1));
    auto T2 = (T_W0[2]) * fpm[2].getTransform(theta_data[i].row(2));

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

#ifdef DEBUG_MODE
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
  std::string prefix;
  if (argc >= 2)
  {
    prefix = argv[1];
    data_input = prefix + "_input_data.txt";
    data_iter = prefix + "_iteration_info.txt";
    data_offset = prefix + "_offset_data.txt";
  }
  else
  {
    data_input = current_workspace + "data/CLOSED_CALIBRATION/1_input_data.txt";
    data_iter = current_workspace + "data/CLOSED_CALIBRATION/1_iteration_info.txt";
    data_offset = current_workspace + "data/CLOSED_CALIBRATION/1__offset_data.txt";    
  }
  std::cout<<"opening: "<<data_input<<std::endl;
  initialize();
  double lambda = 1.0;
  if (argc >= 3)
  {
    std::cout << "\n<<USING SAVED OFFSET>>" << std::endl;
    std::ifstream rf; 
    rf.open(data_offset);
    for (int arm_=0; arm_<N_ARM; arm_++)
    {
      for (int j=0; j<N_J; j++)
      {
        for (int d=0; d<N_DH; d++)
        {
          rf >> offset_matrix[arm_](j,d);
        }
      }
    }
    rf.close();
    std::cout << "offset_matrix LEFT:\n" << offset_matrix[0] << std::endl;
    std::cout << "\noffset_matrix RIGHT:\n" << offset_matrix[1] << std::endl;
    std::cout << "\noffset_matrix TOP:\n" << offset_matrix[2] << "\n\n" <<std::endl;
    lambda = strtod(argv[2],NULL);
    std::cout << "lambda: " << lambda << std::endl;
    data_iter = prefix + "_" + std::to_string(lambda) + "_iteration_info.txt";
    data_offset = prefix + "_" + std::to_string(lambda) + "_offset_data.txt";
  }

  std::ofstream iter_save(data_iter);
  int iter = 100;
  double ev_b;
  double ev_f = 0.0;
  while (iter--)
  {
    std::cout<<"\n\niteration: "<<100-iter<<std::endl;
    for (int i=0; i<N_ARM; i++) {fpm[i].initModel(offset_matrix[i]);}
    ev_b = ev_f;
    p_total = getDistanceDiff();
    ev_f = p_total.squaredNorm() / num_data;
    getJacobian();
    if (iter < 99) {std::cout << "\nrate: " << ((ev_b - ev_f) / ev_b)*100.0 << std::endl;}
    std::cout << "\neval: " << ev_f << std::endl;

    Eigen::Matrix<double, 84, 84> weight;
    weight.setIdentity();
    weight(25,25) = 1e-5;
    weight(25+28*1,25+28*1) = 1e-5;
    weight(25+28*2,25+28*2) = 1e-5;

    auto & j = jacobian;
    
    Eigen::MatrixXd j_diag = (j.transpose() * j ).diagonal().asDiagonal();
    auto j_inv = (j.transpose() * j + lambda * j_diag).inverse() * j.transpose();
    del_phi = weight * j_inv * p_total;
    // std::cout << j.cols() << " row: " << j.rows() << std::endl;
    // std::cout << j << std::endl;
    // std::cout << j_inv.cols() << " row: " << j_inv.rows() << std::endl;
    // std::cout << j_inv << std::endl;
    // std::cout << j_diag.diagonal().transpose() << std::endl;
    // std::cout << j_diag.cols() << " row: " << j_diag.rows() << std::endl;
    // std::cout << j_diag << std::endl;
    // auto j_inv = (j.transpose() * j + lambda * j_diag).inverse() * j.transpose();

    // del_phi = jacobian.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p_total);
    // del_phi = weight * del_phi;

    std::cout << "del_phi.norm(): " << del_phi.norm() << std::endl;
    for (int arm_=0; arm_<N_ARM; arm_++)
    {
      for (int j=0; j<N_J; j++)
      {
        offset_matrix[arm_].row(j) -= del_phi.segment<N_DH>(arm_*N_JDH + j*N_DH);
      }
    }
    std::cout << "\noffset_matrix LEFT:\n" << offset_matrix[0] << std::endl;
    std::cout << "\noffset_matrix RIGHT:\n" << offset_matrix[1] << std::endl;
    std::cout << "\noffset_matrix TOP:\n" << offset_matrix[2] << std::endl;
    if (del_phi.norm() < 1e-9)
    {
      std::cout<<"reached optimal value at iter: "<<100 - iter<<std::endl;
      break;
    }

    iter_save << "iteration: "<< 100 - iter << std::endl;
    iter_save << "eval(before): "<< p_total.squaredNorm() / num_data << std::endl;
    iter_save << "del_phi.norm(): "<< del_phi.norm() << std::endl;
    iter_save << "PANDA LEFT"<< std::endl;
    iter_save << offset_matrix[0].format(tab_format) << std::endl;
    iter_save << "PANDA RIGHT"<< std::endl;
    iter_save << offset_matrix[1].format(tab_format) << std::endl;
    iter_save << "PANDA TOP"<< std::endl;
    iter_save << offset_matrix[2].format(tab_format) << std::endl;
    iter_save <<"\n"<< std::endl;

    std::ofstream offset_save(data_offset);
    offset_save << offset_matrix[0].format(tab_format) <<std::endl;
    offset_save << offset_matrix[1].format(tab_format) <<std::endl;
    offset_save << offset_matrix[2].format(tab_format);
    offset_save.close();
  }

  return 0;
}