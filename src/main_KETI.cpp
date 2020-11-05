#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <mutex>
#include <thread>
#include <robot_model/franka_panda_model.h>
#include "suhan_benchmark.h"

// #define TEST

const int N_DH = 4; // number of dh parameters to calibrate
const int N_J = 7; // number of joints to calibrate, 7 + 2(base+marker)
const int N_CAL = N_DH*N_J; // total variables to calibrate
const int N_XYZ = 3;
const int N_data_input = N_J + N_XYZ;
const std::string arm_ = "panda_top";

std::mutex data_mutex;
std::mutex print_mutex;

double min_cost = 1e100;
Eigen::Matrix<double, N_J, 1> best_sol;
Eigen::Matrix<double, N_CAL, 1> best_dh_sol;

using Vector7d = Eigen::Matrix<double, 7, 1>;
using caldata = std::pair < Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 3, 1> >  ;
// std::vector<Eigen::Matrix<double, N_data_input, 1> > data_input;
std::vector<std::pair< Vector7d, Eigen::Vector3d> > data_input;
Eigen::VectorXd dh_al(N_J), dh_a(N_J), dh_d(N_J);

Eigen::Isometry3d transformKim(const double a, const double d, const double alpha, const double theta)
{
  Eigen::Isometry3d transform_dh;
  transform_dh.setIdentity();
  transform_dh.linear()       << cos(theta),             -1*sin(theta),          0.0,
                                 sin(theta)*cos(alpha),  cos(theta)*cos(alpha),  -1*sin(alpha),
                                 sin(theta)*sin(alpha),  cos(theta)*sin(alpha),  cos(alpha);
  transform_dh.translation()  << a, -1*sin(alpha)*d, cos(alpha)*d;
  return transform_dh;
}

int main(int argc, char**argv)
{
//   srand(time(NULL));
  std::ifstream rf; 
  std::cout << "read data \n"; 
  rf.open(arm_+".txt");
  std::vector <caldata> calib_dataset;
  while (!rf.eof())
  {
    Vector7d q;
    Eigen::Vector3d x;
    // Eigen::Matrix<double, N_data_input, 1> d;
    for (int i=0; i<7; i++)
    {
      rf >> q(i);
    }
    for (int i=0; i<3; i++)
    {
      rf >> x(i);
    }
    q = q / 180.0 * M_PI;
    calib_dataset.push_back(std::make_pair(q,x));
  }
  rf.close();
  std::cout << "complete - size: " << data_input.size() << std::endl;

  if (N_J == 7)
  {
    dh_al << 0.0, -1.0*M_PI_2, M_PI_2, M_PI_2, -1.0*M_PI_2, M_PI_2, M_PI_2;
    dh_a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
    dh_d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;
  }

  else if (N_J == 9)
  {
    dh_al << 0.0, 0.0, -1.0 * M_PI_2, M_PI_2, M_PI_2, -1.0 * M_PI_2, M_PI_2, M_PI_2, 0.0;
    dh_a << 0.0, 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088, 0.0;
    dh_d << 0.0, 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0, 0.107;
  }


  FrankaPandaModel fpm = FrankaPandaModel();

  // base x y z r p y, tip x y z r p y
  auto function_base_tip = [&calib_dataset,&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::VectorXd> out) {
    const auto & q = c.first;
    Eigen::Matrix<double, 7, 4> dh;
    dh.setZero();
    for (int i=0 ; i<7; i++)
    {
      dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);
    }
    Eigen::Matrix<double, 6, 1> x_base, x_tip;
    x_base = x.segment<6>(28);
    x_tip = x.segment<6>(28+6);

    Eigen::Isometry3d t_wb, t_7t;
    t_wb.setIdentity(); t_7t.setIdentity();

    t_wb.translation() = x_base.head<3>();
    t_7t.translation() = x_tip.head<3>();
    t_wb.linear() = (Eigen::AngleAxisd(x_base(3), Eigen::Vector3d::UnitZ()) * 
                     Eigen::AngleAxisd(x_base(4), Eigen::Vector3d::UnitY()) * 
                     Eigen::AngleAxisd(x_base(5), Eigen::Vector3d::UnitX())).matrix();
    t_7t.linear() = (Eigen::AngleAxisd(x_tip(3), Eigen::Vector3d::UnitZ()) * 
                     Eigen::AngleAxisd(x_tip(4), Eigen::Vector3d::UnitY()) * 
                     Eigen::AngleAxisd(x_tip(5), Eigen::Vector3d::UnitX())).matrix();

    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    
    fpm.initModel(dh);
    T_0e = fpm.getTransform(q);

    auto t = (t_wb * T_0e * t_7t).translation();

    out = c.second - t;
  };

  auto jacobian_dh_base_tip = [function_base_tip](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::MatrixXd> out) {
    const auto & q = c.first;
    Eigen::VectorXd y1 = x;
    Eigen::VectorXd y2 = x;
    Eigen::VectorXd t1(N_XYZ);
    Eigen::VectorXd t2(N_XYZ);

    // Use a 7-point central difference stencil on each column.
    for (std::size_t j = 0; j < 28 + 12; j++)
    {
      const double ax = std::fabs(x[j]);
      // Make step size as small as possible while still giving usable accuracy.
      const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);

      // Can't assume y1[j]-y2[j] == 2*h because of precision errors.
      y1[j] += h;
      y2[j] -= h;
      function_base_tip(y1, c, t1);
      function_base_tip(y2, c, t2);
      const Eigen::VectorXd m1 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function_base_tip(y1, c, t1);
      function_base_tip(y2, c, t2);
      const Eigen::VectorXd m2 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function_base_tip(y1, c, t1);
      function_base_tip(y2, c, t2);
      const Eigen::VectorXd m3 = (t1 - t2) / (y1[j] - y2[j]);

      out.col(j) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

      // Reset for next iteration.
      y1[j] = y2[j] = x[j];
    }
  };

  auto function_kim2 = [&calib_dataset,&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::VectorXd> out) {
    const auto & q = c.first;
    Eigen::Matrix<double, N_J, 4> dh;
    dh.setZero();
    for (int i=0 ; i<N_J; i++){dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);}

    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    if (N_J == 7)
    {
      fpm.initModel(dh);
      T_0e = fpm.getTransform(q);
    }
    else if (N_J == 9)
    {
      fpm.initModel9(dh);
      T_0e = fpm.getTransform9(q);
    }

    auto t = T_0e.translation();

    out = c.second - t;
  };

  auto jacobian_dh = [function_kim2](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::MatrixXd> out) {
    const auto & q = c.first;
    Eigen::VectorXd y1 = x;
    Eigen::VectorXd y2 = x;
    Eigen::VectorXd t1(N_XYZ);
    Eigen::VectorXd t2(N_XYZ);

    // Use a 7-point central difference stencil on each column.
    for (std::size_t j = 0; j < N_CAL; j++)
    {
      const double ax = std::fabs(x[j]);
      // Make step size as small as possible while still giving usable accuracy.
      const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);

      // Can't assume y1[j]-y2[j] == 2*h because of precision errors.
      y1[j] += h;
      y2[j] -= h;
      function_kim2(y1, c, t1);
      function_kim2(y2, c, t2);
      const Eigen::VectorXd m1 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function_kim2(y1, c, t1);
      function_kim2(y2, c, t2);
      const Eigen::VectorXd m2 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function_kim2(y1, c, t1);
      function_kim2(y2, c, t2);
      const Eigen::VectorXd m3 = (t1 - t2) / (y1[j] - y2[j]);

      out.col(j) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

      // Reset for next iteration.
      y1[j] = y2[j] = x[j];
    }
  };

  auto jacobian_kim2 = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::MatrixXd> out)
  {
    // const auto & q = c.first;
    Eigen::Matrix<double, 9 ,1> q;
    q.setZero();
    q.segment<7>(1) = c.first;
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie, p_ie_1;
    Eigen::Matrix3d R_0i, R_0i_1;
    Eigen::Isometry3d T_0i, T_0e, T_ie;
    Eigen::Matrix<double, N_XYZ, N_CAL> jacob_k; jacob_k.setZero();
    x_0i << 1.0, 0.0, 0.0; y_0i << 0.0, 1.0, 0.0; z_0i << 0.0, 0.0, 1.0;
    R_0i.setIdentity(); T_0e.setIdentity(); T_0i.setIdentity();
    Eigen::Matrix<double, N_J, 4> dh;
    dh.setZero();
    for (int i=0 ; i<N_J; i++){dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);}
    Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
    Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);
    if (N_J == 7)
    {
      fpm.initModel(dh);
      T_0e = fpm.getTransform(q);
    }
    else if (N_J == 9)
    {
      fpm.initModel9(dh);
      T_0e = fpm.getTransform9(q);
    }

#ifdef TEST
    Eigen::Isometry3d T_0e_test;
    T_0e_test.setIdentity();
    for (int i=0; i<N_J; i++)
    {
      T_0e_test = T_0e_test * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i) + alpha_offset(i), q(i) + q_offset(i));
    }
    if (N_J == 7)
      T_0e_test = T_0e_test * transformKim(0.0, 0.107, 0.0, 0.0);
    else if (N_J == 9)
      T_0e_test = T_0e_test * transformKim(0.0, 0.0, 0.0, 0.0);
    std::cout << "franka model updater:\n" << T_0e.matrix() << std::endl;
    std::cout << "function made:\n" << T_0e_test.matrix() << std::endl;
    std::cout << "translation diff:\n" << (T_0e.translation() - T_0e_test.translation()).norm() << std::endl;
#endif

    p_ie = T_0e.translation();

    for (int i=0; i<N_J; i++)
    { 
      x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i; p_ie_1 = p_ie; R_0i_1 = R_0i;

      x_0i = cos(q(i)+q_offset(i)) * x_0i_1               + sin(q(i)+q_offset(i)) * cos(dh_al(i)+alpha_offset(i)) * y_0i_1 + sin(q(i)+q_offset(i)) * sin(dh_al(i)+alpha_offset(i)) * z_0i_1;
      y_0i = -1.0*sin(q(i)+q_offset(i)) * x_0i_1          + cos(q(i)+q_offset(i)) * cos(dh_al(i)+alpha_offset(i)) * y_0i_1 + cos(q(i)+q_offset(i)) * sin(dh_al(i)+alpha_offset(i)) * z_0i_1;
      z_0i = -1.0*sin(dh_al(i)+alpha_offset(i)) * y_0i_1  + cos(dh_al(i)+alpha_offset(i)) * z_0i_1;

      T_0i = T_0i * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i), dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
      R_0i = T_0i.linear();
      p_ie = (T_0i.inverse() * T_0e).translation();

      jacob_k.col(0 + i*N_DH) += x_0i_1;
      jacob_k.col(1 + i*N_DH) += z_0i;
      jacob_k.col(2 + i*N_DH) += z_0i.cross(R_0i*p_ie);
      if (N_DH == 4){jacob_k.col(3 + i*N_DH) += x_0i_1.cross(R_0i_1*p_ie_1);}
    }
    out = -jacob_k;
  };


  auto print_x = [](const Eigen::Ref<const Eigen::VectorXd> & x) 
  {
    for (int i=0; i<7; i++)
    {
      std::cout << x.segment<4>(i*4).transpose() << std::endl;
    }
    std::cout << x.segment<6>(28).transpose() << std::endl;
    std::cout << x.segment<6>(34).transpose() << std::endl;
  };
  // Eigen::VectorXd x(N_CAL);
  Eigen::VectorXd x(40);
  x.setZero();
#ifdef TEST
  // x = Eigen::Matrix<double, N_CAL, 1>::Random() * 0.01;
  Eigen::MatrixXd j1(N_XYZ,N_CAL), j2(N_XYZ,N_CAL);
  jacobian_dh(x,calib_dataset[1],j1);
  jacobian_kim2(x,calib_dataset[1],j2);
  std::cout << "j1: \n" << j1<<std::endl;
  std::cout << "j2: \n" << j2 << "\ndiff: \n"<< j1-j2 << std::endl;
#else
  Eigen::VectorXd p_total;
  Eigen::MatrixXd jac_total;
  Eigen::VectorXd del_phi;
  int total_len = calib_dataset.size();
  // p_total.resize(N_XYZ*total_len);
  // jac_total.resize(N_XYZ*total_len, N_CAL);
  // del_phi.resize(N_CAL);
  p_total.resize(3*total_len);
  jac_total.resize(3*total_len, 40);
  del_phi.resize(40);
  int iter = 100;
  while (iter--)
  {
    for (int i=0; i<calib_dataset.size(); ++i)
    {
      function_base_tip(x, calib_dataset[i], p_total.segment<3>(i*3));
      std::cout << "p_total[" << i << "]: " << p_total.segment<3>(i*3).transpose() << std::endl;
      std::cout << "calib_dataset[" << i << "]: " << calib_dataset[i].first.transpose() << std::endl;
      std::cout << "calib_dataset[" << i << "]: " << calib_dataset[i].second.transpose() << std::endl;
      jacobian_dh_base_tip(x, calib_dataset[i], jac_total.block<3,40>(i*N_XYZ, 0));
    }
    std::cout << "p_total: \n" << p_total << std::endl;
    std::cout << "J: \n" << jac_total << std::endl;
    std::cout << "\neval: " << p_total.squaredNorm() / total_len << std::endl;
    del_phi = jac_total.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p_total);
    std::cout << "\ndphi:" << std::endl;
    print_x(del_phi);
    x -= del_phi; // jacobi is oppisite direction
    std::cout << "\nx: " << std::endl;
    print_x(x);
    if (del_phi.norm() < 1e-9) break;
  }
  std::ofstream x_out(arm_+"_out.txt");

  Eigen::Matrix<double, N_J, 4> dh;
  dh.setZero();
  for (int i=0 ; i<N_J; i++)
  {
    dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);
  }
  Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
  Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
  Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
  Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);
  // for (auto & q : data_input)
  // {    
  //   Eigen::Isometry3d T_0e;
  //   T_0e.setIdentity();
  //   for (int i=0; i<N_J; i++)
  //   {
  //     T_0e = T_0e * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
  //   }
  //   if (N_J == 7)
  //     T_0e = T_0e * transformKim(0.0, 0.107, 0.0, 0.0);

  //   else if (N_J == 9)
  //     T_0e = T_0e * transformKim(0.0, 0.0, 0.0, 0.0);

  //   auto t = T_0e.translation();
  //   x_out << t.transpose() << std::endl;
  // }
#endif
  return 0;
}