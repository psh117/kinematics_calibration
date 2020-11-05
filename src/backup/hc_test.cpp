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

#define TEST_PRINT

const int N_DH = 4; // number of dh parameters to calibrate
const int N_J = 7; // number of joints to calibrate
const int N_CAL = N_DH*N_J; // total variables to calibrate

struct state_
{
  Eigen::VectorXd dh_alpha(N_J);
  Eigen::VectorXd dh_a(N_J);
  Eigen::VectorXd dh_d(N_J);
  Eigen::VectorXd a_offset;
  Eigen::VectorXd d_offset;
  Eigen::VectorXd theta_offset;
  Eigen::VectorXd alpha_offset;
  Eigen::VectorXd pose_total;
  Eigen::VectorXd del_pose_total;
  Eigen::VectorXd del_phi;
  Eigen::VectorXd EE(3);
  Eigen::VectorXd offset_total(N_CAL); // x(a,d,theta,alpha,a,d,theta,alpha)
  Eigen::Matrix<double, N_J, 4> offset_matrix; // each col (a, d, theta, alpha)
  Eigen::MatrixXd jacobian_total;
  Eigen::Isometry3d T_W0;

  std::vector<Eigen::Matrix<double, N_J, 1>> theta;
  std::string arm_name;
  FrankaPandaModel fpm;
  int num_data;
};

std::map<std::string, std::shared_ptr<a_state_>> arm_state;
state_ left_state_;
state_ right_state_;
state_ top_state_;

Eigen::Isometry3d transformFromDH(const double a, const double d, const double alpha, const double theta)
{
  Eigen::Isometry3d transform_dh;
  transform_dh.setIdentity();
  transform_dh.linear()       << cos(theta),             -1*sin(theta),          0.0,
                                 sin(theta)*cos(alpha),  cos(theta)*cos(alpha),  -1*sin(alpha),
                                 sin(theta)*sin(alpha),  cos(theta)*sin(alpha),  cos(alpha);
  transform_dh.translation()  << a, -1*sin(alpha)*d, cos(alpha)*d;
  return transform_dh;
}

void initTest(state_ &arm)
{
  std::ifstream rf; 
  rf.open(arm.arm_name+".txt");
  while (!rf.eof())
  {
    Eigen::Matrix<double, N_J, 1> d;
    for (int i=0; i<N_J; i++)
    {
      rf >> d(i);
    }
    arm.theta.push_back(d);
  }
  rf.close();
  std::cout << "reading data complete - size: " << arm.theta.size() << std::endl;
  arm.num_data = arm.theta.size();
  arm.pose_total.resize(3*arm.num_data);
  arm.del_pose_total.resize(3*arm.num_data);
  arm.jacobian_total.resize(3*arm.num_data, N_CAL);
  arm.del_phi.resize(N_CAL);

  arm.dh_alpha << 0.0, -1.0*M_PI_2, M_PI_2, M_PI_2, -1.0*M_PI_2, M_PI_2, M_PI_2;
  arm.dh_a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
  arm.dh_d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;
  arm.offset_total.setZero();
  arm.offset_matrix.setZero();
  arm.a_offset.setZero();
  arm.d_offset.setZero();
  arm.theta_offset.setZero();
  arm.alpha_offset.setZero();

  arm.fpm = FrankaPandaModel();
  arm.fpm.initModel(arm.offset_matrix);
}

void updateIter(state_ &arm)
{
  for (int i=0 ; i<N_J; i++){arm.offset_matrix.row(i).head<N_DH>() = arm.offset_total.segment<N_DH>(i*N_DH);}
  arm.a_offset = arm.offset_matrix.col(0);
  arm.d_offset = arm.offset_matrix.col(1);
  arm.theta_offset = arm.offset_matrix.col(2);
  arm.alpha_offset = arm.offset_matrix.col(3);
  arm.fpm.initModel(arm.offset_matrix);
}

void getPose(state_ &arm)
{
  Eigen::Isometry3d T_0e;
  for (int i=0; i<arm.num_data; i++)
  {
    T_0e = fpm.getTransform(arm.theta(i));
    arm.pose_total.segment<3>(i*3) = T_0e.translation();
  }
}

void getJacobian(state_ &arm)
{
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie, p_ie_1;
    Eigen::Matrix3d R_0i, R_0i_1;
    Eigen::Isometry3d T_0i, T_0e, T_ie;
    Eigen::Matrix<double, 3, N_CAL> jacob_k; 

    for (int i=0; i<arm.num_data; i++)
    {
      jacob_k.setZero(); R_0i.setIdentity(); T_0i.setIdentity();
      x_0i << 1.0, 0.0, 0.0; y_0i << 0.0, 1.0, 0.0; z_0i << 0.0, 0.0, 1.0;
      T_0e = fpm.getTransform(arm.theta(i));
      p_ie = T_0e.translation();

      for (int i=0; i<N_J; i++)
      { 
        x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i; p_ie_1 = p_ie; R_0i_1 = R_0i;

        x_0i = cos(arm.theta(i)+arm.theta_offset(i)) * x_0i_1               + sin(arm.theta(i)+arm.theta_offset(i)) * cos(arm.dh_alpha(i)+arm.alpha_offset(i)) * y_0i_1 + sin(arm.theta(i)+arm.theta_offset(i)) * sin(arm.dh_alpha(i)+arm.alpha_offset(i)) * z_0i_1;
        y_0i = -1.0*sin(arm.theta(i)+arm.theta_offset(i)) * x_0i_1          + cos(arm.theta(i)+arm.theta_offset(i)) * cos(arm.dh_alpha(i)+arm.alpha_offset(i)) * y_0i_1 + cos(arm.theta(i)+arm.theta_offset(i)) * sin(arm.dh_alpha(i)+arm.alpha_offset(i)) * z_0i_1;
        z_0i = -1.0*sin(arm.dh_alpha(i)+arm.alpha_offset(i)) * y_0i_1  + cos(arm.dh_alpha(i)+arm.alpha_offset(i)) * z_0i_1;

        T_0i = T_0i * transformFromDH(arm.dh_a(i) + arm.a_offset(i), arm.dh_d(i) + arm.d_offset(i), arm.dh_alpha(i)+arm.alpha_offset(i), arm.theta(i)+arm.theta_offset(i));
        R_0i = T_0i.linear();
        p_ie = (T_0i.inverse() * T_0e).translation();

        jacob_k.col(0 + i*N_DH) += x_0i_1;
        jacob_k.col(1 + i*N_DH) += z_0i;
        jacob_k.col(2 + i*N_DH) += z_0i.cross(R_0i*p_ie);
        if (N_DH == 4){jacob_k.col(3 + i*N_DH) += x_0i_1.cross(R_0i_1*p_ie_1);}
      }

      arm.jacobian_total.block<3,N_CAL>(i*3, 0) = -jacob_k;
    }
}

void getOffset(state_ &arm)
{
    std::cout << "\neval: " << arm.del_pose_total.squaredNorm() / arm.num_data << std::endl;
    arm.del_phi = arm.jacobian_total.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(arm.del_pose_total);
    std::cout << "\ndphi: " << arm.del_phi.transpose() << std::endl;
    arm.offset_total -= arm.del_phi;
    std::cout << "\nx: " << arm.offset_total.transpose() << std::endl;
}

int main(int argc, char**argv)
{
  arm_state.insert(std::pair<std::string,std::shared_ptr<state_>>("panda_left",&left_state_));
  arm_state.insert(std::pair<std::string,std::shared_ptr<state_>>("panda_right",&right_state_));
  arm_state.insert(std::pair<std::string,std::shared_ptr<state_>>("panda_top",&top_state_));

  arm_state["panda_left"] -> EE << 0.0, 0.0, 0.1;
  arm_state["panda_right"] -> EE << 0.0, 0.0, 0.2;
  arm_state["panda_top"] -> EE << 0.0, 0.0, 0.3;

  arm_state["panda_left"] -> T_W0.linear() << 1.0, 0.0, 0.0,
                                              0.0, 1.0, 0.0,
                                              0.0, 0.0, 1.0;
  arm_state["panda_left"] -> T_W0.translation() << 0.0, 0.3, 1.0;

  arm_state["panda_right"] -> T_W0.linear() << 1.0, 0.0, 0.0,
                                               0.0, 1.0, 0.0,
                                               0.0, 0.0, 1.0;
  arm_state["panda_right"] -> T_W0.translation() << 0.0, -0.3, 1.0;

  arm_state["panda_top"] -> T_W0.linear() << -1.0,  0.0, 0.0,
                                              0.0, -1.0, 0.0,
                                              0.0,  0.0, 1.0;
  arm_state["panda_top"] -> T_W0.translation() << 1.35, 0.3, 1.0;

  initTest(arm_state["panda_left"]);
  initTest(arm_state["panda_right"]);
  initTest(arm_state["panda_top"]);

  int iter = 100;
  while (iter--)
  {
    updateIter(*arm_state["panda_left"]);
    getPose(*arm_state["panda_left"]);
    getJacobian(*arm_state["panda_left"]);

    updateIter(*arm_state["panda_right"]);
    getPose(*arm_state["panda_right"]);
    getJacobian(*arm_state["panda_right"]);

    updateIter(*arm_state["panda_top"]);
    getPose(*arm_state["panda_top"]);
    getJacobian(*arm_state["panda_top"]);

    // L -> R -> T -> L
    arm_state["panda_left"] -> del_pose_total = arm_state["panda_right"] -> pose_total - arm_state["panda_left"] -> pose_total;
    arm_state["panda_right"] -> del_pose_total = arm_state["panda_top"] -> pose_total - arm_state["panda_right"] -> pose_total;
    arm_state["panda_top"] -> del_pose_total = arm_state["panda_left"] -> pose_total - arm_state["panda_top"] -> pose_total;

    getOffset(*arm_state["panda_left"])
    getOffset(*arm_state["panda_right"])
    getOffset(*arm_state["panda_top"])

    if (arm_state["panda_left"] -> del_phi.norm() < 1e-9 && arm_state["panda_right"] -> del_phi.norm() < 1e-9 && arm_state["panda_top"] -> del_phi.norm() < 1e-9)
    {
      break;
    }
  }

  updateIter(*arm_state["panda_left"]);
  updateIter(*arm_state["panda_right"]);
  updateIter(*arm_state["panda_top"]);

  std::ofstream left_out("left_out.txt"), right_out("right_out.txt"), top_out("top_out.txt");

  for (auto & q : arm_state["panda_left"] -> theta)
  {    
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<N_J; i++)
    {
      T_0e = T_0e * transformFromDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformFromDH(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_1 << t.transpose() << std::endl;
  }
  return 0;
}