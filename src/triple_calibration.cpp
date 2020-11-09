#include "triple_calibration.h"
// #define DEBUG_MODE

void initialize(const std::string prefix_="1", const bool load_offset_=false)
{

  data_input = current_workspace + prefix_ + "_input_data.txt";
  data_iter = current_workspace + prefix_ + "_iteration_info.txt";
  data_offset = current_workspace + prefix_ + "_offset_data.txt";
  std::cout<<"opening: "<<data_input<<std::endl;

  std::ifstream rf; 
  rf.open(data_input);
  int row = 0;
  while (!rf.eof())
  {
    Eigen::Matrix<double, N_J, 1> d;
    for (auto& arm : arm_state)
    {
      for (int j=0; j<N_J; j++)
      {
        rf >> d(j);
      }
      arm.second -> theta.push_back(d);
    }
    row++;
  }
  rf.close();
  num_data = arm_state["panda_left"] -> theta.size();

  weight.setIdentity();
  weight(25,25) = 1e-5;

  original_dh_matrix << 0.0,     0.333, 0.0,  0.0,
                        0.0,     0.0,   0.0,  -M_PI_2,
                        0.0,     0.316, 0.0,  M_PI_2,
                        0.0825,  0.0,   0.0,  M_PI_2,
                        -0.0825, 0.384, 0.0,  -M_PI_2,
                        0.0,     0.0,   0.0,  M_PI_2,
                        0.088,   0.0,   0.0,  M_PI_2;

  std::cout << "number of datasets: "<<num_data<<std::endl;
  for (auto& arm : arm_state)
  {
    arm.second -> arm_name = arm.first;
    arm.second -> del_pose.resize(N_PARAM*num_data);
    arm.second -> del_pose.setZero();
    arm.second -> del_phi.resize(N_JDH);
    arm.second -> del_phi.setZero();
    arm.second -> jacobian.resize(N_PARAM*num_data, N_JDH);
    arm.second -> jacobian.setZero();
    arm.second -> fpm = FrankaPandaModel();
    arm.second -> offset_matrix.setZero();
    arm.second -> T_EE.setIdentity();
    arm.second -> T_W0.setIdentity();
    if (arm.first == "panda_left")
    {
      arm.second -> arm_true = "panda_right";
      arm.second -> T_W0.translation() << 0.0, 0.3, 1.0;
      arm.second -> T_EE.translation() << 0.1035, -0.1035, 0.0;
    }
    else if (arm.first == "panda_right")
    {
      arm.second -> arm_true = "panda_top";
      arm.second -> T_W0.translation() << 0.0, -0.3, 1.0;
      arm.second -> T_EE.translation() << 0.1035, 0.1035, 0.0;
    }
    else if (arm.first == "panda_top")
    {
      arm.second -> arm_true = "panda_left";
      arm.second -> T_W0.linear() << -1.0,  0.0, 0.0,
                                      0.0, -1.0, 0.0,
                                      0.0,  0.0, 1.0;
      arm.second -> T_W0.translation() << 1.35, 0.3, 1.0;
      arm.second -> T_EE.linear() << -1.0,  0.0, 0.0,
                                      0.0, -1.0, 0.0,
                                      0.0,  0.0, 1.0;
      arm.second -> T_EE.translation() << 0.1035, 0.1035, 0.0;
    }
  }

  if (load_offset_)
  {
    std::cout << "\n<<USING SAVED OFFSET>>" << std::endl;
    std::ifstream rf; 
    rf.open(data_offset);
    for (auto& arm : arm_state)
    {
      for (int j=0; j<N_J; j++)
      {
        for (int d=0; d<N_DH; d++)
        {
          rf >> arm.second -> offset_matrix(j,d);
        }
      }
    std::cout << "offset_matrix " << arm.first << ":\n" << arm.second -> offset_matrix << "\n" << std::endl;
    }
    rf.close();
  }
}

Eigen::Isometry3d transformFromDH(const Eigen::Matrix<double,4,1> dh_)
{
  //// dh_(0) a | dh_(1) d | dh_(2) theta | dh_(3) alpha
  Eigen::Isometry3d transform_dh;
  transform_dh.setIdentity();
  transform_dh.linear()       << cos(dh_(2)),             -1*sin(dh_(2)),          0.0,
                                 sin(dh_(2))*cos(dh_(3)),  cos(dh_(2))*cos(dh_(3)),  -1*sin(dh_(3)),
                                 sin(dh_(2))*sin(dh_(3)),  cos(dh_(2))*sin(dh_(3)),  cos(dh_(3));
  transform_dh.translation()  << dh_(0), -1*sin(dh_(3))*dh_(1), cos(dh_(3))*dh_(1);
  return transform_dh;
}

void getPose(state_ &arm)
{
  for (int i=0; i<num_data; i++)
  {
    Eigen::Isometry3d T_0e = arm.fpm.getTransform(arm.theta[i]);
    Eigen::Isometry3d T_WE = arm.T_W0 * T_0e * arm.T_EE;
    Eigen::Quaterniond quat(T_WE.linear());
    arm.translation.push_back(T_WE.translation());
    arm.quaternion.push_back(quat);
  }
}

void getJacobian(state_ &arm)
{
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie, p_ie_1;
    Eigen::Matrix3d R_0i, R_0i_1;
    Eigen::Isometry3d T_0i, T_0e, T_ie;
    Eigen::Matrix<double, N_PARAM, N_JDH> jacob_k; 

    for (int i=0; i<num_data; i++)
    {
      jacob_k.setZero(); R_0i.setIdentity(); T_0i.setIdentity();
      T_0e = arm.fpm.getTransform(arm.theta[i]);
      p_ie = T_0e.translation();

      for (int i=0; i<N_J; i++)
      { 
        x_0i_1 = R_0i.col(0);
        y_0i_1 = R_0i.col(1);
        z_0i_1 = R_0i.col(2);
        p_ie_1 = T_0i.translation();
        R_0i_1 = R_0i;

        T_0i = T_0i * transformFromDH( original_dh_matrix.row(i) + arm.offset_matrix.row(i) );
        R_0i = T_0i.linear();
        p_ie = (T_0i.inverse() * T_0e).translation();

        jacob_k.block<3,1>(0,0+i*N_DH) = x_0i_1;
        jacob_k.block<3,1>(0,1+i*N_DH) = z_0i;
        jacob_k.block<3,1>(0,2+i*N_DH) = z_0i.cross(R_0i*p_ie);
        if (N_DH == 4){jacob_k.block<3,1>(0,3+i*N_DH) = x_0i_1.cross(R_0i_1*p_ie_1);}
        if (N_PARAM == 6)
        {
          jacob_k.block<3,1>(3,2+i*N_DH) = z_0i;
          if (N_DH == 4){jacob_k.block<3,1>(3,3+i*N_DH) = x_0i_1;}
        }
      }
      arm.jacobian.block<N_PARAM,N_JDH>(i*N_PARAM, 0) = -jacob_k;
    }
}

void getOffset(state_ &arm, state_ &arm_true)
{
    std::cout << "arm: "<< arm.arm_true << std::endl;
    for (int i=0; i<num_data; i++)
    {
      arm.del_pose.segment<3>(i*N_PARAM) = arm_true.translation[i] - arm.translation[i];
      if (N_PARAM > 3)
      { 
        //////////////// NEED TO FIND SOLUTION!!! //////////////////
        arm.del_pose.segment<3>(3+i*N_PARAM) = arm_true.translation[i] - arm.translation[i];
      }
    }

    auto & j = arm.jacobian;
    Eigen::MatrixXd j_diag = (j.transpose() * j ).diagonal().asDiagonal();
    std::cout << "j_diag:\n" << j_diag <<std::endl;
    auto j_i = j_diag;
    auto j_inv = (j.transpose() * j + lambda * (j_diag + j_i.setIdentity())).inverse() * j.transpose();
    arm.del_phi = weight * j_inv * arm.del_pose;
    // std::cout << "j_inv:\n" << j_inv <<std::endl;

    std::cout << "\neval: " << arm.del_pose.squaredNorm() / num_data << std::endl;
    std::cout << "del_phi.norm(): " << arm.del_phi.norm() <<"\n"<< std::endl;

    for (int j=0; j<N_J; j++)
    {
      arm.offset_matrix.row(j) -= arm.del_phi.segment<N_DH>(j*N_DH);
    }
    std::cout << "offset_matrix " << arm.arm_name << ":\n" << arm.offset_matrix << "\n" << std::endl;

}

int main(int argc, char**argv)
{
  state_ left_state_;
  state_ right_state_;
  state_ top_state_;
  arm_state.insert(std::pair<std::string,std::shared_ptr<state_>>("panda_left",&left_state_));
  arm_state.insert(std::pair<std::string,std::shared_ptr<state_>>("panda_right",&right_state_));
  arm_state.insert(std::pair<std::string,std::shared_ptr<state_>>("panda_top",&top_state_));

  if (argc >= 2)
  {
    std::string prefix = argv[1];
    if (argc >= 3)
      initialize(prefix, true);
    else
      initialize(prefix);
  }
  else
    initialize();
  
  std::ofstream iter_save(data_iter);

  int iter = 100;
  while (iter--)
  {
    std::cout<<"iteration: "<<100-iter<<std::endl;
    for (auto& arm : arm_state)
    {
      arm.second -> fpm.initModel(arm.second -> offset_matrix);
      arm.second -> translation.clear();
      arm.second -> quaternion.clear();
      getPose(*arm.second);
      getJacobian(*arm.second);
    }
    for (auto& arm : arm_state) {getOffset(*arm.second, *arm_state[arm.second->arm_true]);}

    if ( arm_state["panda_left"] -> del_phi.norm() + arm_state["panda_right"] -> del_phi.norm() + arm_state["panda_top"] -> del_phi.norm() < 1e-8)
    {
      std::cout<<"reached optimal value at iter: "<<100 - iter<<std::endl;
      break;
    }

    iter_save << "iteration: "<< 100 - iter << std::endl;
    std::ofstream offset_save(data_offset);
    for (auto& arm : arm_state)
    {
      iter_save << "\n" << "arm: "<< arm.first << std::endl;
      iter_save << "eval(before): "<< 100 << std::endl;
      iter_save << "del_phi.norm(): "<< 100 << std::endl;
      if (arm.first == "panda_top")
      {
        offset_save << arm.second -> offset_matrix.format(tab_format);
        iter_save << arm.second -> offset_matrix.format(tab_format) << "\n" << std::endl;
      }
      else
      {
        offset_save << arm.second -> offset_matrix.format(tab_format) << std::endl;
        iter_save << arm.second -> offset_matrix.format(tab_format) << std::endl;
      }
    }
    offset_save.close();
  }

  std::cout<<"testing done."<<std::endl;
  return 0;
}