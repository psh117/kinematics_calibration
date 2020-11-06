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
    arm.second -> pose_total.resize(N_PARAM*num_data);
    arm.second -> pose_total.setZero();
    arm.second -> del_pose_total.resize(N_PARAM*num_data);
    arm.second -> del_pose_total.setZero();
    arm.second -> del_phi.resize(N_PARAM*num_data);
    arm.second -> del_phi.setZero();
    arm.second -> jacobian.resize(N_PARAM*num_data, N_JDH);
    arm.second -> jacobian.setZero();
    arm.second -> fpm = FrankaPandaModel();
    arm.second -> offset_matrix.setZero();
    if (arm.first == "panda_left")
    {
      arm.second -> arm_true = "panda_right";
      arm.second -> T_W0.linear() << 1.0, 0.0, 0.0,
                                     0.0, 1.0, 0.0,
                                     0.0, 0.0, 1.0;
      arm.second -> T_W0.translation() << 0.0, 0.3, 1.0;

    }
    else if (arm.first == "panda_right")
    {
      arm.second -> arm_true = "panda_top";
      arm.second -> T_W0.linear() << 1.0, 0.0, 0.0,
                                     0.0, 1.0, 0.0,
                                     0.0, 0.0, 1.0;
      arm.second -> T_W0.translation() << 0.0, -0.3, 1.0;
    }
    else if (arm.first == "panda_top")
    {
      arm.second -> arm_true = "panda_left";
      arm.second -> T_W0.linear() << -1.0,  0.0, 0.0,
                                      0.0, -1.0, 0.0,
                                      0.0,  0.0, 1.0;
      arm.second -> T_W0.translation() << 1.35, 0.3, 1.0;
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

int main(int argc, char**argv)
{
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

  std::cout<<"testing done."<<std::endl;

  int iter = 100;
  while (iter--)
  {
    std::cout<<"iteration: "<<100-iter<<std::endl;
    for (auto& arm : arm_state)
    {
      arm.second -> fpm.initModel(arm.second -> offset_matrix);
      // calc p_total
      // calc jacobian
      // print eval
      // calc del_phi by LM method
      
    }

    // p_total = getDistanceDiff();
    // getJacobian();

    // std::cout << "\neval: " << p_total.squaredNorm() / num_data << std::endl;

    // auto & j = jacobian;
    // double lambda = 0.01;
    // Eigen::MatrixXd j_diag = (j.transpose() * j ).diagonal().asDiagonal();
    // auto j_inv = (j.transpose() * j + lambda * j_diag).inverse() * j.transpose();
    // del_phi = weight * j_inv * p_total;

    // std::cout << "del_phi.norm(): " << del_phi.norm() <<"\n"<< std::endl;
    // for (int arm_=0; arm_<N_ARM; arm_++)
    // {
    //   for (int j=0; j<N_J; j++)
    //   {
    //     offset_matrix[arm_].row(j) -= del_phi.segment<N_DH>(arm_*N_JDH + j*N_DH);
    //   }
    // }
    // std::cout << "\noffset_matrix LEFT:\n" << offset_matrix[0] << std::endl;
    // std::cout << "\noffset_matrix RIGHT:\n" << offset_matrix[1] << std::endl;
    // std::cout << "\noffset_matrix TOP:\n" << offset_matrix[2] << std::endl;
    // if (del_phi.norm() < 1e-9)
    // {
    //   std::cout<<"reached optimal value at iter: "<<100 - iter<<std::endl;
    //   break;
    // }


    iter_save << "iteration: "<< 100 - iter << std::endl;
    std::ofstream offset_save(data_offset);
    for (auto& arm : arm_state)
    {
      iter_save << "eval(before): "<< 100 << std::endl;
      iter_save << "del_phi.norm(): "<< 100 << std::endl;
      iter_save << arm.first << std::endl;
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
  arm_state.erase("panda_left");
  arm_state.erase("panda_right");
  arm_state.erase("panda_top");
  arm_state.clear();
  return 0;
}