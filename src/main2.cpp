#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <mutex>
#include <thread>
#include <robot_model/franka_panda_model.h>
#include "suhan_benchmark.h"

// // dh test
// int main()
// {
//   auto fpm = FrankaPandaModel();
//   Eigen::VectorXd q(7);
//   q << 0, 0, 0, -1.57, 0, 1.57, 0.77;

//   auto t = fpm.getTransform(q);
//   std::cout << t.matrix() << std::endl;

//   Eigen::Matrix<double, 7,4> dh;
//   dh.setZero();
//   dh(0,0) = 0.01;
//   fpm.initModel(dh);
  

//   auto t2 = fpm.getTransform(q);
//   std::cout << t2.matrix() << std::endl;
//   return 0;
// }

std::mutex data_mutex;
std::mutex print_mutex;

double min_cost = 1e100;
Eigen::Matrix<double, 7, 1> best_sol;
Eigen::Matrix<double, 21, 1> best_dh_sol;

std::vector<Eigen::Matrix<double, 7, 1> > q_input_1, q_input_2, q_input_3;
Eigen::VectorXd dh_al(7), dh_a(7), dh_d(7);
Eigen::Vector3d true_p_1;
Eigen::Vector3d true_p_2;
Eigen::Vector3d true_p_3;

Eigen::Isometry3d transformDH(const double a, const double d, const double alpha, const double theta)
{
  Eigen::Isometry3d transform_dh;
  transform_dh.setIdentity();
  transform_dh.linear()       << cos(theta),             -1*sin(theta),          0.0,
                                 sin(theta)*cos(alpha),  cos(theta)*cos(alpha),  -1*sin(alpha),
                                 sin(theta)*sin(alpha),  cos(theta)*sin(alpha),  cos(alpha);
  transform_dh.translation()  << a, -1*sin(alpha)*d, cos(alpha)*d;

  // Eigen::Isometry3d transform_dh;
  // transform_dh.setIdentity();
  // Eigen::Isometry3d transform_a,transform_d,transform_alpha,transform_theta;
  // transform_a.setIdentity();
  // transform_d.setIdentity();
  // transform_alpha.setIdentity();
  // transform_theta.setIdentity();
  // transform_a.translation() << a, 0.0, 0.0;
  // transform_d.translation() << 0.0, 0.0, d;
  // transform_alpha.linear() << 1.0, 0.0, 0.0, 0.0, cos(alpha), -1.0*sin(alpha), 0.0, sin(alpha), cos(alpha);
  // transform_theta.linear() << cos(theta), -1.0*sin(theta), 0.0, sin(theta), cos(theta), 0.0, 0.0, 0.0, 1.0;
  // transform_dh = transform_alpha * transform_a * transform_d * transform_theta;

  return transform_dh;
}

// void work()
// {
//   auto fpm = FrankaPandaModel();

//   auto function_kim2 = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const Eigen::Ref<const Eigen::VectorXd> &q, Eigen::Ref<Eigen::VectorXd> out) {
//     err.setZero();
//     Eigen::Matrix<double, 7, 4> dh;
//     for (int i=0 ; i<7; i++){dh.row(i).head<3>() = x.segment<3>(i*3);}
//     fpm.initModel(dh);
//     Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
//     auto t = fpm.getTranslation(q+q_offset);
//     out = true_p_1 - t;
//   };

//   auto function_kim = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::VectorXd> out) {
//     Eigen::Vector3d err;
//     err.setZero();
//     Eigen::Matrix<double, 7, 4> dh;
//     for (int i=0 ; i<7; i++)
//     {
//       dh.row(i).head<3>() = x.segment<3>(i*3);
//     }
//     fpm.initModel(dh);
//     Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
//     for (auto & q : q_input_1)
//     {
//       auto t = fpm.getTranslation(q+q_offset);
//       err += true_p_1 - t;
//     }

//     for (auto & q : q_input_2)
//     {
//       auto t = fpm.getTranslation(q+q_offset);
//       err += true_p_2 - t;
//     }
//     for (auto & q : q_input_3)
//     {
//       auto t = fpm.getTranslation(q+q_offset);
//       err += true_p_3 - t;
//     }
//     err /= (q_input_1.size() + q_input_2.size() + q_input_3.size());
//     out = err;
//   };

//   auto jacobian = [function](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out) {
//     Eigen::VectorXd y1 = x;
//     Eigen::VectorXd y2 = x;
//     Eigen::VectorXd t1(1);
//     Eigen::VectorXd t2(1);

//     // Use a 7-point central difference stencil on each column.
//     for (std::size_t j = 0; j < 7; j++)
//     {
//       const double ax = std::fabs(x[j]);
//       // Make step size as small as possible while still giving usable accuracy.
//       const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);

//       // Can't assume y1[j]-y2[j] == 2*h because of precision errors.
//       y1[j] += h;
//       y2[j] -= h;
//       function(y1, t1);
//       function(y2, t2);
//       const Eigen::VectorXd m1 = (t1 - t2) / (y1[j] - y2[j]);
//       y1[j] += h;
//       y2[j] -= h;
//       function(y1, t1);
//       function(y2, t2);
//       const Eigen::VectorXd m2 = (t1 - t2) / (y1[j] - y2[j]);
//       y1[j] += h;
//       y2[j] -= h;
//       function(y1, t1);
//       function(y2, t2);
//       const Eigen::VectorXd m3 = (t1 - t2) / (y1[j] - y2[j]);

//       out.col(j) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

//       // Reset for next iteration.
//       y1[j] = y2[j] = x[j];
//     }
//   };

//   auto jacobian_dh = [function_dh](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out) {
//     Eigen::VectorXd y1 = x;
//     Eigen::VectorXd y2 = x;
//     Eigen::VectorXd t1(1);
//     Eigen::VectorXd t2(1);

//     // Use a 7-point central difference stencil on each column.
//     for (std::size_t j = 0; j < 21; j++)
//     {
//       const double ax = std::fabs(x[j]);
//       // Make step size as small as possible while still giving usable accuracy.
//       const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);

//       // Can't assume y1[j]-y2[j] == 2*h because of precision errors.
//       y1[j] += h;
//       y2[j] -= h;
//       function_dh(y1, t1);
//       function_dh(y2, t2);
//       const Eigen::VectorXd m1 = (t1 - t2) / (y1[j] - y2[j]);
//       y1[j] += h;
//       y2[j] -= h;
//       function_dh(y1, t1);
//       function_dh(y2, t2);
//       const Eigen::VectorXd m2 = (t1 - t2) / (y1[j] - y2[j]);
//       y1[j] += h;
//       y2[j] -= h;
//       function_dh(y1, t1);
//       function_dh(y2, t2);
//       const Eigen::VectorXd m3 = (t1 - t2) / (y1[j] - y2[j]);

//       out.col(j) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

//       // Reset for next iteration.
//       y1[j] = y2[j] = x[j];
//     }
//   };

//   auto jacobian_kim2 = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, const Eigen::Ref<const Eigen::VectorXd> &q, Eigen::Ref<Eigen::MatrixXd> out)
//   {
//     Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie;
//     Eigen::Matrix3d R_0i;
//     Eigen::Isometry3d T_0i, T_0e, T_ie;
//     Eigen::Matrix<double, 3, 21> jacob_k;
//     Eigen::VectorXd alpha(7), a(7), d(7);
//     alpha << 0.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;
//     a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
//     d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;
//     x_0i << 1.0, 0.0, 0.0;
//     y_0i << 0.0, 1.0, 0.0;
//     z_0i << 0.0, 0.0, 1.0;
//     R_0i.setIdentity();
//     Eigen::Matrix<double, 4, 7> dh;
//     for (int i=0 ; i<7; i++){dh.row(i).head<3>() = x.segment<3>(i*3);}
//     fpm.initModel(dh);
//     Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
//     T_0e = fpm.getTransform(q+x);
//     for (int i=0; i<7; i++)
//     {
//       x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
//       if(alpha(i)==0)
//       {
//         x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//         y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//         z_0i = z_0i_1;
//       }
//       else
//       {
//         x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//         y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//         z_0i = alpha(i) * y_0i_1;
//       }
//       T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
//       R_0i = T_0i.linear();
//       p_ie = (T_0i.inverse() * T_0e).translation();
//       jacob_k.col(0 + i*3) += x_0i_1;
//       jacob_k.col(1 + i*3) += z_0i;
//       jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
//       }
//     out = jacob_k;
//     // function_kim2(x, q, d_pos); need (3X1) output
//     // x += (out.transpose() * out).inverse() * out.transpose() * d_pos;
//   }

//   auto jacobian_kim = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out)
//   {
//     Eigen::VectorXd alpha(7), a(7), d(7);
//     alpha << 0.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;
//     a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
//     d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;
//     Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie;
//     Eigen::Matrix3d R_0i;
//     Eigen::Isometry3d T_0i, T_0e, T_ie;
//     Eigen::Matrix<double, 3, 21> jacob_k;
//     x_0i << 1.0, 0.0, 0.0;
//     y_0i << 0.0, 1.0, 0.0;
//     z_0i << 0.0, 0.0, 1.0;
//     R_0i.setIdentity();
//     Eigen::Matrix<double, 4, 7> dh;
//     for (int i=0 ; i<7; i++){dh.row(i).head<3>() = x.segment<3>(i*3);}
//     fpm.initModel(dh);
//     Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
//     for (auto & q : q_input_1)
//     {
//       T_0e = fpm.getTransform(q+x);
//       for (int i=0; i<7; i++)
//       {
//         x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
//         if(alpha(i)==0)
//         {
//           x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//           y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//           z_0i = z_0i_1;
//         }
//         else
//         {
//           x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//           y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//           z_0i = alpha(i) * y_0i_1;
//         }
//         T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
//         R_0i = T_0i.linear();
//         p_ie = (T_0i.inverse() * T_0e).translation();
//         jacob_k.col(0 + i*3) += x_0i_1;
//         jacob_k.col(1 + i*3) += z_0i;
//         jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
//       }
//       T_0i.Identity();
//       x_0i << 1.0, 0.0, 0.0;
//       y_0i << 0.0, 1.0, 0.0;
//       z_0i << 0.0, 0.0, 1.0;
//     }
//     for (auto & q : q_input_2)
//     {
//       T_0e = fpm.getTransform(q+x);
//       for (int i=0; i<7; i++)
//       {
//         x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
//         if(alpha(i)==0)
//         {
//           x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//           y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//           z_0i = z_0i_1;
//         }
//         else
//         {
//           x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//           y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//           z_0i = alpha(i) * y_0i_1;
//         }
//         T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
//         R_0i = T_0i.linear();
//         p_ie = (T_0i.inverse() * T_0e).translation();
//         jacob_k.col(0 + i*3) += x_0i_1;
//         jacob_k.col(1 + i*3) += z_0i;
//         jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
//       }
//       T_0i.Identity();
//       x_0i << 1.0, 0.0, 0.0;
//       y_0i << 0.0, 1.0, 0.0;
//       z_0i << 0.0, 0.0, 1.0;
//     }
//     for (auto & q : q_input_3)
//     {
//       T_0e = fpm.getTransform(q+x);
//       for (int i=0; i<7; i++)
//       {
//         x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
//         if(alpha(i)==0)
//         {
//           x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//           y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//           z_0i = z_0i_1;
//         }
//         else
//         {
//           x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
//           y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
//           z_0i = alpha(i) * y_0i_1;
//         }
//         T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
//         R_0i = T_0i.linear();
//         p_ie = (T_0i.inverse() * T_0e).translation();
//         jacob_k.col(0 + i*3) += x_0i_1;
//         jacob_k.col(1 + i*3) += z_0i;
//         jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
//       }
//       T_0i.Identity();
//       x_0i << 1.0, 0.0, 0.0;
//       y_0i << 0.0, 1.0, 0.0;
//       z_0i << 0.0, 0.0, 1.0;
//     }
//     jacob_k /= (q_input_1.size() + q_input_2.size() + q_input_3.size());
//     out = jacob_k;
//     // function_kim(x, d_pos); need (3X1) output
//     // x += (out.transpose() * out).inverse() * out.transpose() * d_pos;
//   };

//   auto project_dh = [function_dh, jacobian_dh](Eigen::Ref<Eigen::VectorXd> x) {
//     // Newton's method
//     unsigned int iter = 0;
//     double norm = 0;
//     Eigen::VectorXd f(1);
//     Eigen::MatrixXd j(1, 21);

//     Eigen::VectorXd dh_best(21);
//     double best_value = 1e100;

//     function_dh(x, f);
//     int i = 0;
//     double ratio = 1.0;
//     double gamma = 0.95;
//     while (f(0) > 1e-6 && iter++ < 3) 
//     {
//       jacobian_dh(x, j);
//       x -= ratio * j.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(f);
//       ratio *= gamma;
//       function_dh(x, f);
//       data_mutex.lock();
//       if (f(0) < min_cost)
//       {
//         min_cost = f(0);
//         best_dh_sol = x;
//       }
//       data_mutex.unlock();
//       // if (f(0) < best_value)
//       // {
//       //   best_value = f(0);
//       //   q_best = x;
//       // }
//       print_mutex.lock();
//       std::cout << "iter[" << i ++ << "]: " <<  x.transpose() << std::endl
//                 << "    - f: " << f(0) << std::endl;
//       print_mutex.unlock();
//     }
//     x = dh_best;
//     return best_value > 5e-6;
//   };

//   Eigen::VectorXd q(21);
//   q = Eigen::Matrix<double, 7, 1>::Random() * 0.01;
//   Eigen::VectorXd r(1);
//   function_dh(q,r);
//   std::cout << "f: " << r << std::endl;

//   project_dh(q);
//   function_dh(q,r);

// }
int main(int argc, char**argv)
{
//   srand(time(NULL));
  std::ifstream rf; 
  std::cout << "read data 1\n"; 
  rf.open("q_data_1.txt");
  while (!rf.eof())
  {
    Eigen::Matrix<double, 7, 1> d;
    for (int i=0; i<7; i++)
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
    Eigen::Matrix<double, 7, 1> d;
    for (int i=0; i<7; i++)
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
    Eigen::Matrix<double, 7, 1> d;
    for (int i=0; i<7; i++)
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
    Eigen::Matrix<double, 7, 4> dh;
    for (int i=0 ; i<7; i++){dh.row(i).head<3>() = x.segment<3>(i*3);}
    Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
    Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    // Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);
    Eigen::VectorXd alpha_offset(7); alpha_offset.setZero();

    // fpm.initModel(dh);
    // auto t = fpm.getTranslation(q+q_offset);

    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<7; i++){T_0e = T_0e * transformDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));}
    T_0e = T_0e * transformDH(0.0, 0.107, 0.0, 0.0);
    auto t = T_0e.translation();

    out = c.second - t;
  };

  auto jacobian_dh = [function_kim2](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::MatrixXd> out) {
    const auto & q = c.first;
    Eigen::VectorXd y1 = x;
    Eigen::VectorXd y2 = x;
    Eigen::VectorXd t1(3);
    Eigen::VectorXd t2(3);

    // Use a 7-point central difference stencil on each column.
    for (std::size_t j = 0; j < 21; j++)
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
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie;
    Eigen::Matrix3d R_0i;
    Eigen::Isometry3d T_0i, T_0e, T_ie, T_0e_test;
    Eigen::Matrix<double, 3, 21> jacob_k; jacob_k.setZero();
    x_0i << 1.0, 0.0, 0.0; y_0i << 0.0, 1.0, 0.0; z_0i << 0.0, 0.0, 1.0;
    R_0i.setIdentity(); T_0e_test.setIdentity(); T_0i.setIdentity();
    Eigen::Matrix<double, 7, 4> dh;
    for (int i=0 ; i<7; i++){dh.row(i).head<3>() = x.segment<3>(i*3);}
    Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
    Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    // Eigen::Ref<Eigen::VectorXd> alpha_offset = dh.col(3);
    Eigen::VectorXd alpha_offset(7); alpha_offset.setZero();

    for (int i=0; i<7; i++){T_0e_test = T_0e_test * transformDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));}
    T_0e_test = T_0e_test * transformDH(0.0, 0.107, 0.0, 0.0);

    for (int i=0; i<7; i++)
    {
      x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;

      x_0i = cos(q(i)+q_offset(i)) * x_0i_1               + sin(q(i)+q_offset(i)) * cos(dh_al(i)+alpha_offset(i)) * y_0i_1 + sin(q(i)+q_offset(i)) * sin(dh_al(i)+alpha_offset(i)) * z_0i_1;
      y_0i = -1.0*sin(q(i)+q_offset(i)) * x_0i_1          + cos(q(i)+q_offset(i)) * cos(dh_al(i)+alpha_offset(i)) * y_0i_1 + cos(q(i)+q_offset(i)) * sin(dh_al(i)+alpha_offset(i)) * z_0i_1;
      z_0i = -1.0*sin(dh_al(i)+alpha_offset(i)) * y_0i_1  + cos(dh_al(i)+alpha_offset(i)) * z_0i_1;

      T_0i = T_0i * transformDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i), dh_al(i), q(i)+q_offset(i));
      R_0i = T_0i.linear();
      p_ie = (T_0i.inverse() * T_0e_test).translation();

      jacob_k.col(0 + i*3) += x_0i_1;
      jacob_k.col(1 + i*3) += z_0i;
      jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
    }
    out = jacob_k;
    // function_kim2(x, q, d_pos); need (3X1) output
    // x += (out.transpose() * out).inverse() * out.transpose() * d_pos;
  };

  Eigen::VectorXd q(7);
  Eigen::VectorXd x(21);
  q << 0, 0, 0, -1.57, 0, 1.57, 0.79;
  // q << 0, 0, 0.79, -1.57, 0, 1.57, 0.79;
  // x = Eigen::Matrix<double, 21, 1>::Random() * 0.01;
  x.setZero();
  Eigen::VectorXd r1(1), r2(1);
  Eigen::MatrixXd j1(3,21), j2(3,21);
  // jacobian_dh(x,calib_dataset[0],j1);
  // jacobian_kim2(x,calib_dataset[0],j2);
  // std::cout << "j1: \n" << j1<<std::endl;
  // std::cout << "j2: \n" << j2 << "\ndiff: \n"<< j1+j2 << std::endl;

  Eigen::VectorXd p_total;
  Eigen::MatrixXd jac_total;
  Eigen::VectorXd del_phi;
  int total_len = q_input_1.size() + q_input_2.size() + q_input_3.size();
  p_total.resize(3*total_len);
  jac_total.resize(3*total_len, 21);
  del_phi.resize(21);

  int iter = 100;
  while (iter--)
  {
    for (int i=0; i<calib_dataset.size(); ++i)
    {
      function_kim2(x, calib_dataset[i], p_total.segment<3>(i*3));
      jacobian_kim2(x, calib_dataset[i], jac_total.block<3,21>(i*3, 0));
    }
    std::cout << "eval: " << p_total.squaredNorm() / total_len << std::endl;
    del_phi = jac_total.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(p_total);
    std::cout << "dphi: " << del_phi.transpose() << std::endl;
    x += del_phi; // jacobi is oppisite direction
    std::cout << "x: " << x.transpose() << std::endl;
    if (del_phi.norm() < 1e-9) break;
  }
  std::ofstream x_out_1("x_out_1.txt"), x_out_2("x_out_2.txt"), x_out_3("x_out_3.txt");

  Eigen::Matrix<double, 7, 4> dh;
  for (int i=0 ; i<7; i++)
  {
    dh.row(i).head<3>() = x.segment<3>(i*3);
  }
  Eigen::Ref<Eigen::VectorXd> a_offset = dh.col(0);
  Eigen::Ref<Eigen::VectorXd> d_offset = dh.col(1);
  Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
  Eigen::VectorXd alpha_offset(7); 
  alpha_offset.setZero();
  // auto fpm = FrankaPandaModel();
  for (auto & q : q_input_1)
  {    
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<7; i++)
    {
      T_0e = T_0e * transformDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformDH(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_1 << t.transpose() << std::endl;
  }
  for (auto & q : q_input_2)
  {
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<7; i++)
    {
      T_0e = T_0e * transformDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformDH(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_2 << t.transpose() << std::endl;
  }
  for (auto & q : q_input_3)
  {
    Eigen::Isometry3d T_0e;
    T_0e.setIdentity();
    for (int i=0; i<7; i++)
    {
      T_0e = T_0e * transformDH(dh_a(i) + a_offset(i), dh_d(i) + d_offset(i) ,  dh_al(i)+alpha_offset(i), q(i)+q_offset(i));
    }
    T_0e = T_0e * transformDH(0.0, 0.107, 0.0, 0.0);

    auto t = T_0e.translation();
    x_out_3 << t.transpose() << std::endl;
  }

  return 0;
}