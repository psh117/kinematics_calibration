

  auto function_closed = [&calib_dataset,&fpm_l, &fpm_r](const Eigen::Ref<const Eigen::VectorXd> &x, const caldata &c, Eigen::Ref<Eigen::VectorXd> out) {
    const auto & q = c.first;
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

    Eigen::Isometry3d t_b71, t_b72;
    fpm_l.initModel(dh.topRows<N_SJ>());
    fpm_r.initModel(dh.bottomRows<N_SJ>());
    t_b71 = fpm_l.getTransform(c.head<7>());
    t_b72 = fpm_r.getTransform(c.tail<7>());
    auto t = T_0e.translation();

    out = c.second - t;
  };
