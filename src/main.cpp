#include <iostream>
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

Eigen::Vector3d true_p_1;
Eigen::Vector3d true_p_2;
Eigen::Vector3d true_p_3;

Eigen::Matrix4d transformDH(const double a, const double d, const double alpha, const double theta)
{
  Eigen::Matrix4d dh_x;
  dh_x << 1.0,        0.0,            0.0,            a,
          0.0,        cos(alpha),     -1*sin(alpha),  0.0,
          0.0,        sin(alpha),     cos(alpha),     0.0,
          0.0,        0.0,            0.0,            1.0;
  Eigen::Matrix4d dh_z;
  dh_z << cos(theta), -1*sin(theta),  0.0,            0.0,
          sin(theta), cos(theta),     0.0,            0.0,
          0.0,        0.0,            1.0,            d,
          0.0,        0.0,            0.0,            1.0;
  return dh_x * dh_z;
}

void work()
{
  auto fpm = FrankaPandaModel();

  auto function_dh = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::VectorXd> out) {
    double err = 0.0;
    Eigen::Matrix<double, 7, 4> dh;
    for (int i=0 ; i<7; i++)
    {
      dh.row(i).head<3>() = x.segment<3>(i*3);
    }
    fpm.initModel(dh);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    for (auto & q : q_input_1)
    {
      auto t = fpm.getTranslation(q+q_offset);
      err += (t-true_p_1).squaredNorm();
    }

    for (auto & q : q_input_2)
    {
      auto t = fpm.getTranslation(q+q_offset);
      err += (t-true_p_2).squaredNorm();
    }
    for (auto & q : q_input_3)
    {
      auto t = fpm.getTranslation(q+q_offset);
      err += (t-true_p_3).squaredNorm();
    }
    err /= (q_input_1.size() + q_input_2.size() + q_input_3.size());

    out(0) = err;
  };

  auto function_kim = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::VectorXd> out) {
    Eigen::Vector3d err;
    err.setZero();
    Eigen::Matrix<double, 7, 4> dh;
    for (int i=0 ; i<7; i++)
    {
      dh.row(i).head<3>() = x.segment<3>(i*3);
    }
    fpm.initModel(dh);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    for (auto & q : q_input_1)
    {
      auto t = fpm.getTranslation(q+q_offset);
      err += true_p_1 - t;
    }

    for (auto & q : q_input_2)
    {
      auto t = fpm.getTranslation(q+q_offset);
      err += true_p_2 - t;
    }
    for (auto & q : q_input_3)
    {
      auto t = fpm.getTranslation(q+q_offset);
      err += true_p_3 - t;
    }
    err /= (q_input_1.size() + q_input_2.size() + q_input_3.size());
    out = err;
  };

  auto function = [&fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::VectorXd> out) {
    double err = 0.0;
    for (auto & q : q_input_1)
    {
      auto t = fpm.getTranslation(q+x);
      err += (t-true_p_1).squaredNorm();
    }

    for (auto & q : q_input_2)
    {
      auto t = fpm.getTranslation(q+x);
      err += (t-true_p_2).squaredNorm();
    }
    for (auto & q : q_input_3)
    {
      auto t = fpm.getTranslation(q+x);
      err += (t-true_p_3).squaredNorm();
    }
    err /= (q_input_1.size() + q_input_2.size() + q_input_3.size());

    out(0) = err;
  };
  auto jacobian = [function](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out) {
    Eigen::VectorXd y1 = x;
    Eigen::VectorXd y2 = x;
    Eigen::VectorXd t1(1);
    Eigen::VectorXd t2(1);

    // Use a 7-point central difference stencil on each column.
    for (std::size_t j = 0; j < 7; j++)
    {
      const double ax = std::fabs(x[j]);
      // Make step size as small as possible while still giving usable accuracy.
      const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);

      // Can't assume y1[j]-y2[j] == 2*h because of precision errors.
      y1[j] += h;
      y2[j] -= h;
      function(y1, t1);
      function(y2, t2);
      const Eigen::VectorXd m1 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function(y1, t1);
      function(y2, t2);
      const Eigen::VectorXd m2 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function(y1, t1);
      function(y2, t2);
      const Eigen::VectorXd m3 = (t1 - t2) / (y1[j] - y2[j]);

      out.col(j) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

      // Reset for next iteration.
      y1[j] = y2[j] = x[j];
    }
  };

  auto jacobian_dh = [function_dh](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out) {
    Eigen::VectorXd y1 = x;
    Eigen::VectorXd y2 = x;
    Eigen::VectorXd t1(1);
    Eigen::VectorXd t2(1);

    // Use a 7-point central difference stencil on each column.
    for (std::size_t j = 0; j < 21; j++)
    {
      const double ax = std::fabs(x[j]);
      // Make step size as small as possible while still giving usable accuracy.
      const double h = std::sqrt(std::numeric_limits<double>::epsilon()) * (ax >= 1 ? ax : 1);

      // Can't assume y1[j]-y2[j] == 2*h because of precision errors.
      y1[j] += h;
      y2[j] -= h;
      function_dh(y1, t1);
      function_dh(y2, t2);
      const Eigen::VectorXd m1 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function_dh(y1, t1);
      function_dh(y2, t2);
      const Eigen::VectorXd m2 = (t1 - t2) / (y1[j] - y2[j]);
      y1[j] += h;
      y2[j] -= h;
      function_dh(y1, t1);
      function_dh(y2, t2);
      const Eigen::VectorXd m3 = (t1 - t2) / (y1[j] - y2[j]);

      out.col(j) = 1.5 * m1 - 0.6 * m2 + 0.1 * m3;

      // Reset for next iteration.
      y1[j] = y2[j] = x[j];
    }
  };

  auto jacobian2 = [function, &fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out)
  {
    Eigen::Matrix<double, 1, 7> grad;
    grad.setZero(1,7);
    for (auto & q : q_input_1)
    {
      auto t = fpm.getTranslation(q+x);
      grad += (t-true_p_1).transpose() * fpm.getJacobianMatrix(q+x).block<3,7>(0,0);
    }

    for (auto & q : q_input_2)
    {
      auto t = fpm.getTranslation(q+x);
      grad += (t-true_p_2).transpose() * fpm.getJacobianMatrix(q+x).block<3,7>(0,0);
    }
    for (auto & q : q_input_3)
    {
      auto t = fpm.getTranslation(q+x);
      grad += (t-true_p_3).transpose() * fpm.getJacobianMatrix(q+x).block<3,7>(0,0);
    }

    grad /= (q_input_1.size() + q_input_2.size() + q_input_3.size()) / 2.0;

    out = grad; //function3d(x).transpose() * fpm.getJacobianMatrix(x).block<3,7>(0,0);
  };

  auto jacobian_kim = [function_dh, &fpm](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::MatrixXd> out)
  {
    Eigen::VectorXd d_pos(3);
    Eigen::VectorXd alpha(7), a(7), d(7);
    alpha << 0.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0;
    a << 0.0, 0.0, 0.0, 0.0825, -0.0825, 0.0, 0.088;
    d << 0.333, 0.0, 0.316, 0.0, 0.384, 0.0, 0.0;
    Eigen::Vector3d x_0i_1, y_0i_1, z_0i_1, x_0i, y_0i, z_0i, p_ie;
    Eigen::Matrix3d R_0i;
    Eigen::Isometry3d T_0i, T_0e, T_ie;
    Eigen::Matrix<double, 3, 21> jacob_k;
    x_0i << 1.0, 0.0, 0.0;
    y_0i << 0.0, 1.0, 0.0;
    z_0i << 0.0, 0.0, 1.0;
    R_0i.setIdentity();
    Eigen::Matrix<double, 4, 7> dh;
    for (int i=0 ; i<7; i++)
    {
      dh.row(i).head<3>() = x.segment<3>(i*3);
    }
    fpm.initModel(dh);
    Eigen::Ref<Eigen::VectorXd> q_offset = dh.col(2);
    for (auto & q : q_input_1)
    {
      T_0e = fpm.getTransform(q+x);
      for (int i=0; i<7; i++) // calculate J_a, J_alpha, J_theta
      {
        x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
        if(alpha(i)==0)
        {
          x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
          y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
          z_0i = z_0i_1;
        }
        else
        {
          x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
          y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
          z_0i = alpha(i) * y_0i_1;
        }
        T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
        R_0i = T_0i.linear();
        p_ie = (T_0i.inverse() * T_0e).translation();
        jacob_k.col(0 + i*3) += x_0i_1;
        jacob_k.col(1 + i*3) += z_0i;
        jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
      }
      T_0i.Identity();
      x_0i << 1.0, 0.0, 0.0;
      y_0i << 0.0, 1.0, 0.0;
      z_0i << 0.0, 0.0, 1.0;
    }
    for (auto & q : q_input_2)
    {
      T_0e = fpm.getTransform(q+x);
      for (int i=0; i<7; i++) // calculate J_a, J_alpha, J_theta
      {
        x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
        if(alpha(i)==0)
        {
          x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
          y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
          z_0i = z_0i_1;
        }
        else
        {
          x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
          y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
          z_0i = alpha(i) * y_0i_1;
        }
        T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
        R_0i = T_0i.linear();
        p_ie = (T_0i.inverse() * T_0e).translation();
        jacob_k.col(0 + i*3) += x_0i_1;
        jacob_k.col(1 + i*3) += z_0i;
        jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
      }
      T_0i.Identity();
      x_0i << 1.0, 0.0, 0.0;
      y_0i << 0.0, 1.0, 0.0;
      z_0i << 0.0, 0.0, 1.0;
    }
    for (auto & q : q_input_3)
    {
      T_0e = fpm.getTransform(q+x);
      for (int i=0; i<7; i++) // calculate J_a, J_alpha, J_theta
      {
        x_0i_1 = x_0i; y_0i_1 = y_0i; z_0i_1 = z_0i;
        if(alpha(i)==0)
        {
          x_0i = sin(q(i)+x(2 + i*3)) * y_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
          y_0i = cos(q(i)+x(2 + i*3)) * y_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
          z_0i = z_0i_1;
        }
        else
        {
          x_0i = alpha(i) * sin(q(i)+x(2 + i*3)) * z_0i_1 + cos(q(i)+x(2 + i*3)) * x_0i_1;
          y_0i = alpha(i) * cos(q(i)+x(2 + i*3)) * z_0i_1 + sin(q(i)+x(2 + i*3)) * x_0i_1;
          z_0i = alpha(i) * y_0i_1;
        }
        T_0i = T_0i * transformDH(a(i) + x(0 + i*3), d(i) + x(1 + i*3), M_PI_2 * alpha(i), q(i)+x(2 + i*3));
        R_0i = T_0i.linear();
        p_ie = (T_0i.inverse() * T_0e).translation();
        jacob_k.col(0 + i*3) += x_0i_1;
        jacob_k.col(1 + i*3) += z_0i;
        jacob_k.col(2 + i*3) += z_0i.cross(R_0i*p_ie);
      }
      T_0i.Identity();
      x_0i << 1.0, 0.0, 0.0;
      y_0i << 0.0, 1.0, 0.0;
      z_0i << 0.0, 0.0, 1.0;
    }
    jacob_k /= (q_input_1.size() + q_input_2.size() + q_input_3.size());
    out = jacob_k;
    // function_kim(x, d_pos); need (3X1) output
    // x += (out.transpose() * out).inverse() * out.transpose() * d_pos;
  };

  auto project = [function, jacobian2](Eigen::Ref<Eigen::VectorXd> x) {
    // Newton's method
    unsigned int iter = 0;
    double norm = 0;
    Eigen::VectorXd f(1);
    Eigen::MatrixXd j(1, 7);

    Eigen::VectorXd q_best(7);
    double best_value = 1e100;

    function(x, f);
    int i = 0;
    double ratio = 1.0;
    double gamma = 0.95;
    while (f(0) > 1e-6 && iter++ < 15) 
    {
      jacobian2(x, j);
      x -= ratio * j.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(f);
      ratio *= gamma;
      function(x, f);
      data_mutex.lock();
      if (f(0) < min_cost)
      {
        min_cost = f(0);
        best_sol = x;
      }
      data_mutex.unlock();
      // if (f(0) < best_value)
      // {
      //   best_value = f(0);
      //   q_best = x;
      // }
      print_mutex.lock();
      std::cout << "iter[" << i ++ << "]: " <<  x.transpose() << std::endl
                << "    - f: " << f(0) << std::endl;
      print_mutex.unlock();
    }
    x = q_best;
    return best_value > 5e-6;
  };
  auto project_dh = [function_dh, jacobian_dh](Eigen::Ref<Eigen::VectorXd> x) {
    // Newton's method
    unsigned int iter = 0;
    double norm = 0;
    Eigen::VectorXd f(1);
    Eigen::MatrixXd j(1, 21);

    Eigen::VectorXd dh_best(21);
    double best_value = 1e100;

    function_dh(x, f);
    int i = 0;
    double ratio = 1.0;
    double gamma = 0.95;
    while (f(0) > 1e-6 && iter++ < 3) 
    {
      jacobian_dh(x, j);
      x -= ratio * j.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(f);
      ratio *= gamma;
      function_dh(x, f);
      data_mutex.lock();
      if (f(0) < min_cost)
      {
        min_cost = f(0);
        best_dh_sol = x;
      }
      data_mutex.unlock();
      // if (f(0) < best_value)
      // {
      //   best_value = f(0);
      //   q_best = x;
      // }
      print_mutex.lock();
      std::cout << "iter[" << i ++ << "]: " <<  x.transpose() << std::endl
                << "    - f: " << f(0) << std::endl;
      print_mutex.unlock();
    }
    x = dh_best;
    return best_value > 5e-6;
  };

  Eigen::VectorXd q(21);
  q = Eigen::Matrix<double, 7, 1>::Random() * 0.01;
  Eigen::VectorXd r(1);
  function_dh(q,r);
  std::cout << "f: " << r << std::endl;

  project_dh(q);
  function_dh(q,r);

}
int main(int argc, char**argv)
{
  srand(time(NULL));
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

  auto fpm = FrankaPandaModel();

  std::vector<std::thread> processes;

  for (int i = 0; i < 8; i++)
  {
    processes.push_back(std::thread(work));
  }

  for (auto & process : processes)
  {
    process.join();
  }

  std::cout << "f: " << min_cost << std::endl;

  std::cout << best_sol.transpose() << std::endl;
  std::cout << "in deg: " << best_sol.transpose() * 180.0 / 3.141592653589 << std::endl;

  std::ofstream x_out_1("x_out_1.txt"), x_out_2("x_out_2.txt"), x_out_3("x_out_3.txt");

  // auto fpm = FrankaPandaModel();
  for (auto & q_raw : q_input_1)
  {
    auto t = fpm.getTranslation(q_raw + best_sol);
    x_out_1 << t.transpose() << std::endl;
  }
  for (auto & q_raw : q_input_2)
  {
    auto t = fpm.getTranslation(q_raw + best_sol);
    x_out_2 << t.transpose() << std::endl;
  }
  for (auto & q_raw : q_input_3)
  {
    auto t = fpm.getTranslation(q_raw + best_sol);
    x_out_3 << t.transpose() << std::endl;
  }

  return 0;
}