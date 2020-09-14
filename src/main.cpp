#include <iostream>
#include <fstream>
#include <vector>
#include <robot_model/franka_panda_model.h>

int main()
{
  auto fpm = FrankaPandaModel();
  std::ifstream rf; 
  std::vector<Eigen::Matrix<double, 7, 1> > q_input_1, q_input_2;
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

  Eigen::Vector3d true_p_1;
  Eigen::Vector3d true_p_2;

  true_p_1 << -0.55, 0.0, 0.068;
  true_p_2 << 0.0, -0.45, 0.068;

  auto function = [&fpm, &q_input_1, &q_input_2, &true_p_1, &true_p_2](const Eigen::Ref<const Eigen::VectorXd> &x, Eigen::Ref<Eigen::VectorXd> out) {
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
    err /= (q_input_1.size() + q_input_2.size());

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

  auto project = [function, jacobian](Eigen::Ref<Eigen::VectorXd> x) {
    // Newton's method
    unsigned int iter = 0;
    double norm = 0;
    Eigen::VectorXd f(1);
    Eigen::MatrixXd j(1, 7);

    Eigen::VectorXd q_best(7);
    double best_value = 1e100;

    function(x, f);
    int i = 0;
    double ratio = 1;
    double gamma = 0.95;
    while (f(0) > 1e-6 && iter++ < 50) 
    {
      jacobian(x, j);
      x -= ratio * j.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(f);
      ratio *= gamma;
      function(x, f);
      if (f(0) < best_value)
      {
        best_value = f(0);
        q_best = x;
      }
      std::cout << "iter[" << i ++ << "]: " <<  x.transpose() << std::endl
                << "    - f: " << f(0) << std::endl;
    }
    x = q_best;
    return best_value > 5e-6;
  };

  Eigen::VectorXd q(7);
  q.setZero();
  Eigen::VectorXd r(1);
  function(q,r);
  std::cout << "f: " << r << std::endl;

  project(q);
  function(q,r);
  std::cout << "f: " << r << std::endl;

  std::cout << q.transpose() << std::endl;
  std::cout << "in deg: " << q.transpose() * 180.0 / 3.141592653589 << std::endl;

  std::ofstream x_out_1("x_out_1.txt"), x_out_2("x_out_2.txt");

  for (auto & q_raw : q_input_1)
  {
    auto t = fpm.getTranslation(q_raw + q);
    x_out_1 << t.transpose() << std::endl;
  }
  for (auto & q_raw : q_input_2)
  {
    auto t = fpm.getTranslation(q_raw + q);
    x_out_2 << t.transpose() << std::endl;
  }

  return 0;
}