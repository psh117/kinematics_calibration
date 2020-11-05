#pragma once

#include <memory>

#include <rbdl/rbdl.h>
#include <Eigen/Dense>

#include "robot_model/robot_model.h"


class FrankaPandaModel : public RobotModel
{
public:
  static constexpr int kDof=7;
  static constexpr int newDof=9;
  static constexpr double kGravity=9.8;

  FrankaPandaModel();

  Eigen::MatrixXd getJacobianMatrix(const Eigen::VectorXd &q) const override;
  Eigen::Vector3d getTranslation(const Eigen::VectorXd &q) const override;
  Eigen::Matrix3d getRotation(const Eigen::VectorXd &q) const override;
  Eigen::Isometry3d getTransform(const Eigen::VectorXd &q) const override;
  Eigen::MatrixXd getJointLimit() const override;
  Eigen::VectorXd getInitialConfiguration() const override;
  Eigen::Vector3d getTranslation9(const Eigen::VectorXd &q) const;
  Eigen::Matrix3d getRotation9(const Eigen::VectorXd &q) const;
  Eigen::Isometry3d getTransform9(const Eigen::VectorXd &q) const;

  int getDof() override;
  void initModel();
  void initModel(const Eigen::Ref<const Eigen::Matrix<double, 7, 4> > & dh); // 7x4
  void initModel9(const Eigen::Ref<const Eigen::Matrix<double, 9, 4> > & dh); // 9x4

private:
  std::shared_ptr<RigidBodyDynamics::Model> model_;
  RigidBodyDynamics::Math::Vector3d ee_position_;
  Eigen::Matrix3d rot_ee_;

  RigidBodyDynamics::Math::Vector3d com_position_[kDof];
  RigidBodyDynamics::Joint joint_[kDof];
  Eigen::Vector3d joint_posision_[kDof];
  RigidBodyDynamics::Body body_[kDof];
  unsigned int body_id_[kDof];

  RigidBodyDynamics::Math::Vector3d com_position9_[newDof];
  RigidBodyDynamics::Joint joint9_[newDof];
  Eigen::Vector3d joint_posision9_[newDof];
  RigidBodyDynamics::Body body9_[newDof];
  unsigned int body_id9_[newDof];
};
