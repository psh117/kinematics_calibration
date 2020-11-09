#define TEST_PRINT

void getPose(state_ &arm)
{
  Eigen::Isometry3d T_0e;
  for (int i=0; i<arm.num_data; i++)
  {
    T_0e = fpm.getTransform(arm.theta(i));
    arm.pose_total.segment<3>(i*3) = T_0e.translation();
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

  return 0;
}