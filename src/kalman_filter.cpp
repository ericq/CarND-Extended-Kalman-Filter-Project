#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
 	x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  // x_ reprents predicted state after the above Predict() invocation. 
  // X_ contains px,py,vx,vy. H vector muliplification is just a fancy way to 
  // abandons vx,vy. thus z_pred only contains px, py
 	VectorXd z_pred = H_ * x_; 
 	VectorXd y = z - z_pred;
 	MatrixXd Ht = H_.transpose();
 	MatrixXd S = H_ * P_ * Ht + R_;
 	MatrixXd Si = S.inverse();
 	MatrixXd PHt = P_ * Ht;
 	MatrixXd K = PHt * Si;

 	//new estimate
 	x_ = x_ + (K * y);
 	long x_size = x_.size();
 	MatrixXd I = MatrixXd::Identity(x_size, x_size);
 	P_ = (I - K * H_) * P_;
}

/**
 *  input is a vector contains px, py, vx, vy
 */
VectorXd KalmanFilter::CartesianToPolar(const VectorXd &cartesianInput) {
  VectorXd polarOutput (3);
  double px = cartesianInput(0), py = cartesianInput(1), vx = cartesianInput(2), vy = cartesianInput(3);
  
  polarOutput << sqrt(px*px + py*py), atan(py/px), (px*vx+py*vy)/sqrt(px*px + py*py);

  return polarOutput;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  // Jacobian Hj to repace H_, compared to the standard KF
  Tools tool;
  
  MatrixXd Hj = tool.CalculateJacobian(x_);

  VectorXd z_pred = Hj * x_;

  // convert z from Cartesian coordindates to Polar coordindates
  //VectorXd h_z(3);
  //h_z << CartesianToPolar(z_pred);

 	VectorXd y = z - z_pred;
 	MatrixXd Hjt = Hj.transpose();
 	MatrixXd S = Hj * P_ * Hjt + R_;
 	MatrixXd Si = S.inverse();
 	MatrixXd PHjt = P_ * Hjt;
 	MatrixXd K = PHjt * Si;

 	//new estimate
 	x_ = x_ + (K * y);
 	long x_size = x_.size();
 	MatrixXd I = MatrixXd::Identity(x_size, x_size);
 	P_ = (I - K * Hj) * P_;
}
