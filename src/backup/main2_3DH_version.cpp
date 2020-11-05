#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <mutex>
#include <thread>
#include <robot_model/franka_panda_model.h>
#include "suhan_benchmark.h"

const int N_DH = 3; // number of dh parameters to calibrate
const int N_J = 7; // number of joints to calibrate
const int N_CAL = N_DH*N_J; // total variables to calibrate

std::mutex data_mutex;
std::mutex print_mutex;

double min_cost = 1e100;
Eigen::Matrix<double, N_J, 1> best_sol;
Eigen::Matrix<double, N_CAL, 1> best_dh_sol;

std::vector<Eigen::Matrix<double, N_J, 1> > q_input_1, q_input_2, q_input_3;
Eigen::VectorXd dh_al(N_J), dh_a(N_J), dh_d(N_J);
Eigen::Vector3d true_p_1;
Eigen::Vector3d true_p_2;
Eigen::Vector3d true_p_3;

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
  std::cout << "read data 1\n"; 
  rf.open("q_data_1.txt");
  while (!rf.eof())
  {
    Eigen::Matrix<double, N_J, 1> d;
    for (int i=0; i<N_J; i++)
    {
      rf >> d(i);
    }
    q_input_1.push_back(d);
  }
  rf.close();
  std::cout << "complete - size: " << q_input_1.size() << std::endl;

  std::cout << "read data 2\n";
  rf.open("q_data_2.txt");
  while (!rf.eof())
  {
    Eigen::Matrix<double, N_J, 1> d;
    for (int i=0; i<N_J; i++)
    {
      rf >> d(i);
    }
    q_input_2.push_back(d);
  }
  rf.close();
  std::cout << "complete - size: " << q_input_2.size() << std::endl;

  std::cout << "read data 3\n";
  rf.open("q_data_3.txt");
  while (!rf.eof())
  {
    Eigen::Matrix<double, N_J, 1> d;
    for (int i=0; i<N_J; i++)
    {
      rf >> d(i);
    }
    q_input_3.push_back(d);
  }
  rf.close();
  std::cout << "complete - size: " << q_input_3.size() << std::endl;
  
  double z = 0.065;
  if (argc  == 2)
  {
    z = std::atof(argv[1]); 
  }
  true_p_1 << -0.65, 0.0, z;
  true_p_2 << 0.0, 0.45, z;
  true_p_3 << 0.675, 0, z;

  typedef std::pair < Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 3, 1> > caldata ;

  std::vector <caldata> calib_dataset;
  // std::map < Eigen::Matrix<double, 7, 1>, Eigen::Matrix<double, 3, 1>
  // memory prealloc
  calib_dataset.reserve(q_input_1.size() + q_input_2.size() + q_input_3.size());

  for (auto & q : q_input_1)
  {
    calib_dataset.push_back(std::make_pair(q, true_p_1));
  }
  for (auto & q : q_input_2)
  {
    calib_dataset.push_back(std::make_pair(q, true_p_2));
  }
  for (auto & q : q_input_3)
  {
    calib_dataset.push_back(std::make_pair(q, true_p_3));
  }

  dh_al << 0.0, -1.0*M_PI_2, M_PI_2, M_PI_2, -1.0*M_PI_2, M_PI_2, M_PI_2;
  dh_a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
  dh_d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;

  auto fpm = FrankaPandaModel();

  auto function_kim2 = [&calib_dataset,&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::VectorXd> out) {
    const auto & q = c.first;
    Eigen::Matrix<double, N_J, 4> dh;
    dh.setZero();
    for (int i=0 ; i<N_J; i++){dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);}
    Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
    Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    // Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);
    Eigen::VectorXd alpha_offset(7); alpha_offset.setZero();

    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    Eigen::VectorXd alpha(7), a(7), d(7);
    for (int i=0; i<7; i++)
    {
      T_0e = T_0e * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformKim(0.0, 0.107, 0.0, 0.0);
    auto t = T_0e.translation();

    // fpm.initModel(dh);
    // auto t = fpm.getTranslation(q+q_offset);

    out = c.second - t;
  };

  auto jacobian_dh = [function_kim2](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::MatrixXd> out) {
    const auto & q = c.first;
    Eigen::VectorXd y1 = x;
    Eigen::VectorXd y2 = x;
    Eigen::VectorXd t1(3);
    Eigen::VectorXd t2(3);

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
    const auto & q = c.first;
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie, p_ie_1;
    Eigen::Matrix3d R_0i, R_0i_1;
    Eigen::Isometry3d T_0i, T_0e, T_ie;
    Eigen::Matrix<double, 3, N_CAL> jacob_k; jacob_k.setZero();
    x_0i << 1.0, 0.0, 0.0; y_0i << 0.0, 1.0, 0.0; z_0i << 0.0, 0.0, 1.0;
    R_0i.setIdentity(); T_0e.setIdentity(); T_0i.setIdentity();
    Eigen::Matrix<double, N_J, 4> dh;
    dh.setZero();
    for (int i=0 ; i<N_J; i++){dh.row(i).head<N_DH>() = x.segment<N_DH>(i*N_DH);}
    Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
    Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    // Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);
    Eigen::VectorXd alpha_offset(7); alpha_offset.setZero();

    fpm.initModel(dh);
    T_0e = fpm.getTransform(q+q_offset);
    p_ie = T_0e.translation();

    Eigen::Isometry3d T_0e_test;
    T_0e_test.setIdentity();
    for (int i=0; i<7; i++)
    {
      T_0e_test = T_0e_test * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e_test = T_0e_test * transformKim(0.0, 0.107, 0.0, 0.0);

    // std::cout << "franka model updater:\n" << T_0e.matrix() << std::endl;
    // std::cout << "function made:\n" << T_0e_test.matrix() << std::endl;
    // std::cout << "translation diff:\n" << (T_0e.translation() - T_0e_test.translation()).norm() << std::endl;

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

  Eigen::VectorXd q(N_J);
  Eigen::VectorXd x(N_CAL);
  q << 0, 0, 0, -1.57, 0, 1.57, 0.79;
  // q << 0, 0, 0.79, -1.57, 0, 1.57, 0.79;
  // x = Eigen::Matrix<double, N_CAL, 1>::Random() * 0.01;
  x.setZero();
  Eigen::VectorXd r1(1), r2(1);
  Eigen::MatrixXd j1(3,N_CAL), j2(3,N_CAL);
  // jacobian_dh(x,calib_dataset[0],j1);
  // jacobian_kim2(x,calib_dataset[0],j2);
  // std::cout << "j1: \n" << j1<<std::endl;
  // std::cout << "j2: \n" << j2 << "\ndiff: \n"<< j1-j2 << std::endl;

  Eigen::VectorXd p_total;
  Eigen::MatrixXd jac_total;
  Eigen::VectorXd del_phi;
  int total_len = q_input_1.size() + q_input_2.size() + q_input_3.size();
  p_total.resize(3*total_len);
  jac_total.resize(3*total_len, N_CAL);
  del_phi.resize(N_CAL);

  int iter = 100;
  while (iter--)
  {
    for (int i=0; i<calib_dataset.size(); ++i)
    {
      function_kim2(x, calib_dataset[i], p_total.segment<3>(i*3));
      jacobian_kim2(x, calib_dataset[i], jac_total.block<3,N_CAL>(i*3, 0));
    }
    std::cout << "eval: " << p_total.squaredNorm() / total_len << std::endl;
    del_phi = jac_total.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p_total);
    std::cout << "dphi: " << del_phi.transpose() << std::endl;
    x += del_phi; // jacobi is oppisite direction
    std::cout << "x: " << x.transpose() << std::endl;
    if (del_phi.norm() < 1e-9) break;
  }
  std::ofstream x_out_1("x_out_1.txt"), x_out_2("x_out_2.txt"), x_out_3("x_out_3.txt");

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
  // auto fpm = FrankaPandaModel();
  for (auto & q : q_input_1)
  {    
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<N_J; i++)
    {
      T_0e = T_0e * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformKim(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_1 << t.transpose() << std::endl;
  }
  for (auto & q : q_input_2)
  {
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<N_J; i++)
    {
      T_0e = T_0e * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformKim(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_2 << t.transpose() << std::endl;
  }
  for (auto & q : q_input_3)
  {
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<N_J; i++)
    {
      T_0e = T_0e * transformKim(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformKim(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_3 << t.transpose() << std::endl;
  }

  return 0;
}