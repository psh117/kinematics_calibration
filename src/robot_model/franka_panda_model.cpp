#include "robot_model/franka_panda_model.h"
#include <Eigen/Geometry>


FrankaPandaModel::FrankaPandaModel()
{
  initModel();
}

Eigen::MatrixXd FrankaPandaModel::
getJacobianMatrix(const Eigen::VectorXd &q) const
{
  Eigen::MatrixXd j_temp;
  j_temp.setZero(6,kDof);
  Eigen::Matrix<double, 6, kDof> j;
  RigidBodyDynamics::CalcPointJacobian6D(*model_, q, body_id_[kDof - 1], ee_position_, j_temp, true);

  for (int i = 0; i<2; i++)
  {
    j.block<3, kDof>(i * 3, 0) = j_temp.block<3, kDof>(3 - i * 3, 0);
  }

  return j;
}

Eigen::Vector3d FrankaPandaModel::getTranslation(const Eigen::VectorXd &q) const
{
  return RigidBodyDynamics::CalcBodyToBaseCoordinates(*model_, q, body_id_[kDof - 1], ee_position_, true);
}

Eigen::Matrix3d FrankaPandaModel::getRotation(const Eigen::VectorXd &q) const
{
  Eigen::Matrix3d rotation = RigidBodyDynamics::CalcBodyWorldOrientation
      (*model_, q, body_id_[kDof - 1], true).transpose();
  Eigen::Matrix3d body_to_ee_rotation;
  body_to_ee_rotation.setIdentity();
  body_to_ee_rotation(1, 1) = -1;
  body_to_ee_rotation(2, 2) = -1;
  // return Eigen::AngleAxisd(M_PI/4., Eigen::Vector3d::UnitZ()) * body_to_ee_rotation * rotation;
  // return rotation * Eigen::AngleAxisd(M_PI/4., Eigen::Vector3d::UnitZ()) * body_to_ee_rotation;
  return rotation * body_to_ee_rotation;
}

Eigen::Isometry3d FrankaPandaModel::getTransform(const Eigen::VectorXd &q) const
{
  Eigen::Isometry3d transform;
  transform.linear() = getRotation(q);
  transform.translation() = getTranslation(q);

  return transform;
}

Eigen::MatrixXd FrankaPandaModel::getJointLimit() const
{
  Eigen::Matrix<double, kDof, 2> joint_limits;
  joint_limits << -2.8973	,	2.8973	,
                  -1.7628	,	1.7628	,
                  -2.8973	,	2.8973	,
                  -3.0718	,	-0.0698	,
                  -2.8973	,	2.8973	,
                  -0.0175	,	3.7525	,
                  -2.8973	,	2.8973;
  return joint_limits;
}

Eigen::VectorXd FrankaPandaModel::getInitialConfiguration() const
{
}

int FrankaPandaModel::getDof()
{
  return kDof;
}

void FrankaPandaModel::initModel()
{
  model_ = std::make_shared<RigidBodyDynamics::Model>();

  model_->gravity = Eigen::Vector3d(0., 0., -kGravity);

  double mass[kDof];
  mass[0] = 1.0;
  mass[1] = 1.0;
  mass[2] = 1.0;
  mass[3] = 1.0;
  mass[4] = 1.0;
  mass[5] = 1.0;
  mass[6] = 1.0;

  Eigen::Vector3d axis[kDof];
  axis[0] = Eigen::Vector3d::UnitZ();
  axis[1] = Eigen::Vector3d::UnitY();
  axis[2] = Eigen::Vector3d::UnitZ();
  axis[3] = -1.0*Eigen::Vector3d::UnitY();
  axis[4] = Eigen::Vector3d::UnitZ();
  axis[5] = -1.0*Eigen::Vector3d::UnitY();
  axis[6] = -1.0*Eigen::Vector3d::UnitZ();

  Eigen::Vector3d global_joint_position[kDof];

  global_joint_position[0] = Eigen::Vector3d(0.0, 0.0, 0.3330);
  global_joint_position[1] = global_joint_position[0];
  global_joint_position[2] = Eigen::Vector3d(0.0, 0.0, 0.6490);
  global_joint_position[3] = Eigen::Vector3d(0.0825, 0.0, 0.6490);
  global_joint_position[4] = Eigen::Vector3d(0.0, 0.0, 1.0330);
  global_joint_position[5] = Eigen::Vector3d(0.0, 0.0, 1.0330);
  global_joint_position[6] = Eigen::Vector3d(0.0880, 0.0, 1.0330);

  joint_posision_[0] = global_joint_position[0];
  for (int i = 1; i < kDof; i++)
    joint_posision_[i] = global_joint_position[i] - global_joint_position[i - 1];

  com_position_[0] = Eigen::Vector3d(0.000096, -0.0346, 0.2575);
  com_position_[1] = Eigen::Vector3d(0.0002, 0.0344, 0.4094);
  com_position_[2] = Eigen::Vector3d(0.0334, 0.0266, 0.6076);
  com_position_[3] = Eigen::Vector3d(0.0331, -0.0266, 0.6914);
  com_position_[4] = Eigen::Vector3d(0.0013, 0.0423, 0.9243);
  com_position_[5] = Eigen::Vector3d(0.0421, -0.0103, 1.0482);
  com_position_[6] = Eigen::Vector3d(0.1, -0.0120, 0.9536);
  ee_position_ = Eigen::Vector3d(0.0880, 0, 0.9260);
  ee_position_ -= global_joint_position[6];

  for (int i = 0; i < kDof; i++)
    com_position_[i] -= global_joint_position[i];

    RigidBodyDynamics::Math::Vector3d inertia[kDof];
  for (int i = 0; i < kDof; i++)
    inertia[i] = Eigen::Vector3d::Identity() * 0.001;

    for (int i = 0; i < kDof; i++) {
        body_[i] = RigidBodyDynamics::Body(mass[i], com_position_[i], inertia[i]);
        joint_[i] = RigidBodyDynamics::Joint(RigidBodyDynamics::JointTypeRevolute, axis[i]);
        if (i == 0)
            body_id_[i] = model_->AddBody(0, RigidBodyDynamics::Math::Xtrans(joint_posision_[i]), joint_[i], body_[i]);
        else
            body_id_[i] = model_->AddBody(body_id_[i - 1], RigidBodyDynamics::Math::Xtrans(joint_posision_[i]), joint_[i], body_[i]);
    }
}


void FrankaPandaModel::initModel(const Eigen::Ref<const Eigen::Matrix<double, 7, 4> > & dh)
{
  model_ = std::make_shared<RigidBodyDynamics::Model>();

  model_->gravity = Eigen::Vector3d(0., 0., -kGravity);

  double mass[kDof];
  mass[0] = 1.0;
  mass[1] = 1.0;
  mass[2] = 1.0;
  mass[3] = 1.0;
  mass[4] = 1.0;
  mass[5] = 1.0;
  mass[6] = 1.0;

  Eigen::Vector3d axis[kDof];
  axis[0] = Eigen::Vector3d::UnitZ();
  axis[1] = Eigen::Vector3d::UnitY();
  axis[2] = Eigen::Vector3d::UnitZ();
  axis[3] = -1.0*Eigen::Vector3d::UnitY();
  axis[4] = Eigen::Vector3d::UnitZ();
  axis[5] = -1.0*Eigen::Vector3d::UnitY();
  axis[6] = -1.0*Eigen::Vector3d::UnitZ();

  Eigen::Vector3d global_joint_position[kDof];

  // global_joint_position[0] = Eigen::Vector3d(0.0 + dh(0,0), 0.0, 0.3330 + dh(0,1));
  // global_joint_position[1] = global_joint_position[0];
  // global_joint_position[1](2) += dh(1,1);
  // global_joint_position[1](2) += dh(1,1);
  // global_joint_position[2] = Eigen::Vector3d(0.0, 0.0, 0.6490);
  // global_joint_position[3] = Eigen::Vector3d(0.0825, 0.0, 0.6490);
  // global_joint_position[4] = Eigen::Vector3d(0.0, 0.0, 1.0330);
  // global_joint_position[5] = Eigen::Vector3d(0.0, 0.0, 1.0330);
  // global_joint_position[6] = Eigen::Vector3d(0.0880, 0.0, 1.0330);

  // joint_posision_[0] = global_joint_position[0];
  // for (int i = 1; i < kDof; i++)
  //   joint_posision_[i] = global_joint_position[i] - global_joint_position[i - 1];
  
  const Eigen::Ref<const Eigen::VectorXd> &a = dh.col(0);
  const Eigen::Ref<const Eigen::VectorXd> &d = dh.col(1);
  using v3 = Eigen::Vector3d;
  //                      ---- x -----    ---- y ----    ---- z ----
  joint_posision_[0] = v3( 0.0    + a(0), 0.0,        0.3330 + d(0));
  joint_posision_[1] = v3( 0.0    + a(1), 0.0 + d(1), 0.0          );
  joint_posision_[2] = v3( 0.0    + a(2), 0.0,        0.316  + d(2));
  joint_posision_[3] = v3( 0.0825 + a(3), 0.0 - d(3), 0.0          );
  joint_posision_[4] = v3(-0.0825 + a(4), 0.0,        0.384  + d(5));
  joint_posision_[5] = v3( 0.0    + a(5), 0.0 - d(5), 0.0          );
  joint_posision_[6] = v3( 0.088  + a(6), 0.0,        - d(6)        );
  ee_position_ = v3(0.0, 0.0, -0.107);
  com_position_[0] = Eigen::Vector3d(0.000096, -0.0346, 0.2575);
  com_position_[1] = Eigen::Vector3d(0.0002, 0.0344, 0.4094);
  com_position_[2] = Eigen::Vector3d(0.0334, 0.0266, 0.6076);
  com_position_[3] = Eigen::Vector3d(0.0331, -0.0266, 0.6914);
  com_position_[4] = Eigen::Vector3d(0.0013, 0.0423, 0.9243);
  com_position_[5] = Eigen::Vector3d(0.0421, -0.0103, 1.0482);
  com_position_[6] = Eigen::Vector3d(0.1, -0.0120, 0.9536);
  // ee_position_ = Eigen::Vector3d(0.0880, 0, 0.9260);
  // ee_position_ -= global_joint_position[6];

  for (int i = 0; i < kDof; i++)
    com_position_[i] -= global_joint_position[i];

    RigidBodyDynamics::Math::Vector3d inertia[kDof];
  for (int i = 0; i < kDof; i++)
    inertia[i] = Eigen::Vector3d::Identity() * 0.001;

    for (int i = 0; i < kDof; i++) {
        body_[i] = RigidBodyDynamics::Body(mass[i], com_position_[i], inertia[i]);
        joint_[i] = RigidBodyDynamics::Joint(RigidBodyDynamics::JointTypeRevolute, axis[i]);
        if (i == 0)
            body_id_[i] = model_->AddBody(0, RigidBodyDynamics::Math::Xtrans(joint_posision_[i]), joint_[i], body_[i]);
        else
            body_id_[i] = model_->AddBody(body_id_[i - 1], RigidBodyDynamics::Math::Xtrans(joint_posision_[i]), joint_[i], body_[i]);
    }
}

