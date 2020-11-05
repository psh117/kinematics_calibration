
    switch (state_)
    {
      case EXEC :
        run_time = time.toSec() - motion_start_time_;
        force = f_ext.head<3>() - accumulated_wrench_.head<3>();

        f_star = PegInHole::threeDofMove(origin, current_, target_pos, xd, time.toSec(), motion_start_time_, duration_);
        m_star = PegInHole::keepCurrentOrientation(origin, current_, xd, 200, 5);

        if (run_time > 0.05 && Criteria::checkContact(force, Eigen::Isometry3d::Identity(), contact_force_))
        {
          std::cout << "CHECK CONTATCT!!!!!" << std::endl;
          succeed_flag++;
          state_ = KEEPCURRENT;
          is_mode_changed_ = true;
        }
        else if ((origin.translation() - current_.translation()).norm >= limit_dist_)
        {
          std::cout << "DISTANCE REACHED!!!!!" << std::endl;
          succeed_flag++;
          state_ = KEEPCURRENT;
          is_mode_changed_ = true;
        }
        else if (timeOut(time.toSec(), arm.task_start_time_.toSec(), duration_))
        {
          std::cout << "Time out" << std::endl;
          setAborted();
        }
        break;

      case KEEPCURRENT:
        f_star = PegInHole::keepCurrentPose(origin, current_, xd, 5000, 200, 200, 5).head<3>();
        m_star = PegInHole::keepCurrentOrientation(origin, current_, xd, 200, 5);
        break;
    }



    switch (state_)
    {
      case EXEC :
        run_time = time.toSec() - motion_start_time_;
        force = f_ext.head<3>() - accumulated_wrench_.head<3>();

        f_star = PegInHole::threeDofMove(origin, current_, target_pos, xd, time.toSec(), motion_start_time_, duration_);
        m_star = PegInHole::keepCurrentOrientation(origin, current_, xd, 200, 5);

        if (run_time > 0.05 && Criteria::checkContact(force, Eigen::Isometry3d::Identity(), contact_force_))
        {
          std::cout << "CHECK CONTATCT!!!!!" << std::endl;
          succeed_flag++;
          state_ = KEEPCURRENT;
          is_mode_changed_ = true;
        }
        else if ((origin.translation() - current_.translation()).norm >= limit_dist_)
        {
          std::cout << "DISTANCE REACHED!!!!!" << std::endl;
          succeed_flag++;
          state_ = KEEPCURRENT;
          is_mode_changed_ = true;
        }
        else if (timeOut(time.toSec(), arm.task_start_time_.toSec(), duration_))
        {
          std::cout << "Time out" << std::endl;
          setAborted();
        }
        break;

      case KEEPCURRENT:
        f_star = PegInHole::keepCurrentPose(origin, current_, xd, 5000, 200, 200, 5).head<3>();
        m_star = PegInHole::keepCurrentOrientation(origin, current_, xd, 200, 5);
        break;
    }


bool AssembleTripleMoveActionServer::compute(ros::Time time)
{
  if (!control_running_)
    return false;
  if (!as_.isActive())
    return false;
  if (mu_.find("panda_left") != mu_.end() && mu_.find("panda_right") != mu_.end())
  {
    if (first_)
    {
      motion_start_time_ = time.toSec();
      first_ = false;
    }
    computeArm(time, *mu_["panda_left"], left_arm_origin_, left_target, left_state_);
    computeArm(time, *mu_["panda_right"], right_arm_origin_, right_target, right_state_);
    computeArm(time, *mu_["panda_top"], top_arm_origin_, top_target, top_state_);
    if (succeed_flag >= 2)
      setSucceeded();
    return true;
  }
  else
  {
    ROS_ERROR("[AssembleTripleMoveActionServer::goalCallback] the name are not in the arm list.");
    return false;
  }
}