#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  if (!is_initialized_) {
    x_.fill(0.0);
    x_(1) = -0.05;
    n_x_ = 5;
    n_aug_ = 7;
    P_ = MatrixXd::Identity(n_x_,n_x_);
    P_(0,0) = 0.1;
    P_(1,1) = 0.1;

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double phidot = meas_package.raw_measurements_(2);

      x_(0) = rho*cos(phi);
      x_(1) = rho*sin(phi);
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }


  delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
  Prediction(delta_t);
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }

  time_us_ = meas_package.timestamp_;
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /// Augmented state and covariance
  lambda_ = 3 - n_aug_;

  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  VectorXd x_aug_ = VectorXd(n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_ + 1);

  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = pow(std_a_,2.0);
  P_aug(6,6) = pow(std_yawdd_,2.0);

  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug_;

  for(int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_)*A.col(i);
    Xsig_aug.col(i + n_aug_ + 1) = x_aug_ - sqrt(lambda_ + n_aug_)*A.col(i);
  }


  /// Predicted mean and covariance
  for(int i = 0; i < 2*n_aug_ + 1; i++){
    double px = Xsig_aug.col(i)(0);
    double py = Xsig_aug.col(i)(1);
    double v = Xsig_aug.col(i)(2);
    double phi = Xsig_aug.col(i)(3);
    double phidot = Xsig_aug.col(i)(4);
    double nu_a = Xsig_aug.col(i)(5);
    double nu_phidot = Xsig_aug.col(i)(6);

    double px_p,py_p;

    if(fabs(phidot) > 0.001){
      px_p = px + v/phidot*(sin(phi + phidot*delta_t) - sin(phi));
      py_p = py + v/phidot*(cos(phi) - cos(phi + phidot*delta_t));
    }
    else{
      px_p = px + v*delta_t*cos(phi);
      py_p = py + v*delta_t*sin(phi);
    }

    double v_p, phi_p, phidot_p;

    v_p = v;
    phi_p = phi + phidot*delta_t;
    phidot_p  = phidot;


    px_p += 0.5*nu_a*pow(delta_t,2.0)*cos(phi);
    py_p += 0.5*nu_a*pow(delta_t,2.0)*sin(phi);
    v_p += nu_a*delta_t;

    phi_p += 0.5*nu_phidot*pow(delta_t,2.0);
    phidot_p += nu_phidot*delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = phi_p;
    Xsig_pred_(4,i) = phidot_p;
  }

  /// Predict state and covariance matrix
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_/(lambda_ + n_aug_);

  for(int i = 1; i < 2*n_aug_ + 1; i++){
    weights_(i) = 1/(2*(lambda_+n_aug_));
  }
  x_.fill(0.0);
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++){
    x_+= weights_(i)*Xsig_pred_.col(i);
  }

  for(int i = 0; i < 2*n_aug_ + 1; i++){
    VectorXd xdiff(n_x_);
    xdiff = Xsig_pred_.col(i) - x_;
    xdiff(3) = atan2(sin(xdiff(3)),cos(xdiff(3)));
    P_ += weights_(i)*xdiff*xdiff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  MatrixXd H_laser(2,5);
          H_laser << 1, 0, 0, 0, 0,
                     0, 1, 0, 0, 0;

  VectorXd z_pred = H_laser*x_;
  VectorXd z(2);

  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  MatrixXd R_(2,2);
  R_ << pow(std_laspx_,2.0), 0,
          0, pow(std_laspy_,2.0);

  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser.transpose();
  MatrixXd S = H_laser*P_*Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_*Ht;
  MatrixXd K = PHt*Si;

  x_ += K*y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  P_ = (I - K*H_laser)*P_;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /// Predicted measurement
  int n_z = 3; /// Degrees of freedom for radar
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S = MatrixXd(n_z,n_z);

  VectorXd z = VectorXd(n_z);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  for(int i = 0; i < 2*n_aug_ + 1; i++){
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double phi = Xsig_pred_(3,i);
    double phidot = Xsig_pred_(4,i);


    Zsig(0,i) = sqrt(pow(px,2.0) + pow(py,2.0));
    Zsig(1,i) = atan2(py,px);
    Zsig(2,i) = (px*cos(phi)*v + py*sin(phi)*v)/Zsig(0,i);
  }

  z_pred.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++){
    z_pred += weights_(i)*Zsig.col(i);
  }

  S.fill(0.0);

  MatrixXd R = MatrixXd(n_z,n_z);
  R << pow(std_radr_,2.0), 0, 0,
       0, pow(std_radphi_,2.0), 0,
       0, 0, pow(std_radrd_,2.0);

  for(int i = 0; i < 2*n_aug_ + 1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while(z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)< M_PI) z_diff(1)+=2.*M_PI;
    S += weights_(i)*z_diff*z_diff.transpose();
  }

  S += R;

  /// Update state, gain, covariance

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < 2*n_aug_ + 1; i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)-=2.*M_PI;

    Tc += weights_(i)*x_diff*z_diff.transpose();
  }

  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc*S.inverse();

  x_ += K*(z - z_pred);
  P_ -= K*S*K.transpose();
}
