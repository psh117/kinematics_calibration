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

#define DEBUG_MODE

std::string arm_1 = "panda_left";
std::string arm_2 = "panda_right";
std::string arm_3 = "panda_top";
std::string user_name = "kimhc";
std::string current_workspace = "/home/" + user_name + "/git/kinematics_calibration/";
std::string data_input = current_workspace + "data/CLOSED_CALIBRATION/q.txt";
std::string middle_input = current_workspace + "data/CLOSED_CALIBRATION/3_input_data.txt";
std::string data_output = current_workspace + "data/CLOSED_CALIBRATION/collected_q.txt";
Eigen::IOFormat tab_format(Eigen::FullPrecision, 0, "\t", "\n");

const int N_ROBOT = 3;  //number of robots
const int N_J = 7;    // number of joints in one robot
const int N_TOT = N_ROBOT * N_J;    // total joints
const double DEL_MOVE = 1e-3 * 3.0;

void collect_one()
{
  int input_counter_ = 0; int middle_counter_ = 0;
  Eigen::VectorXd q_1(N_TOT); q_1.setZero();
  std::vector<Eigen::VectorXd> saver_;
  std::ifstream rf; 
  rf.open(data_input);
  while (!rf.eof())
  { 
    input_counter_++;
    Eigen::VectorXd new_q_1(N_TOT);

    for(int j=0; j<N_TOT; j++)
    {
      rf >> new_q_1(j);
    }

    if((q_1-new_q_1).norm() > DEL_MOVE)
    {
      middle_counter_++;
      saver_.push_back(new_q_1);
      q_1 = new_q_1;
      
    }
  }
  rf.close();
  std::ofstream middle_work(middle_input);
  for (int i=0; i<saver_.size();i++)
  {
    if (i==saver_.size()-1)
      middle_work << saver_[i].transpose().format(tab_format);
    else
      middle_work << saver_[i].transpose().format(tab_format)<<std::endl;
  }
  middle_work.close();
#ifdef DEBUG_MODE
  std::cout<<"input size: "<<input_counter_<<"\nmiddle size: "<<middle_counter_<<std::endl;
#endif
}

void collect_two()
{
  int collect_counter_ = 0;
  std::vector<Eigen::VectorXd> collector;
  bool is_first_ = true;
  std::ifstream rf; 
  rf.open(middle_input);
  while (!rf.eof())
  { 
    Eigen::VectorXd new_q_1(N_TOT);

    for(int j=0; j<N_TOT; j++)
    {
      rf >> new_q_1(j);
    }

    if (is_first_)
    {
      is_first_ = false;
      collector.push_back(new_q_1);
      collect_counter_++;
    }
    bool is_save_ = true;
    for(int i=0; i<collector.size(); i++)
    {
      if((collector[i]-new_q_1).norm() < DEL_MOVE)
      {
        is_save_ = false;
      }
      if (i == collector.size()-1)
      {
        if (is_save_)
        {
          collect_counter_++;
          collector.push_back(new_q_1);
        }
        is_save_ = true;
      }

    }
  }
  std::ofstream output_work(data_output);
  for (int i=0; i<collector.size();i++)
  {
    if (i == collector.size()-1)
      output_work << collector[i].transpose().format(tab_format);
    else
      output_work << collector[i].transpose().format(tab_format)<<std::endl;
  }
  rf.close();
  output_work.close();
#ifdef DEBUG_MODE
  std::cout<<"collected size: "<<collect_counter_<<std::endl;
#endif
}

int main(int argc, char**argv)
{
#ifdef DEBUG_MODE
  std::cout<<"reading from -- "<<data_input<<std::endl;
#endif
  collect_one();
  // collect_two();

  return 0;
}