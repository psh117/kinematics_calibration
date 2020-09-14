#include <iostream>
#include <robot_model/franka_panda_model.h>

int main()
{
  auto fpm = FrankaPandaModel();

  Eigen::VectorXd q(7);
  q.setZero();
  auto t = fpm.getTransform(q);

  std::cout << t.matrix() << std::endl;
  return 0;
}