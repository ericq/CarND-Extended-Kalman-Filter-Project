#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

 	//create a 4D state vector, we don't know yet the values of the x state
  ekf_.x_ = VectorXd(4);

	//state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;  
  

  //measurement matrix
	ekf_.H_ = MatrixXd(2, 4);
	ekf_.H_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

	//the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;


 	//set the acceleration noise components
	noise_ax = 5;
	noise_ay = 5;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    

  //cout<< "entering ProcessMeasurement()... is_initialized_= " << is_initialized_ << endl;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    //cout << "entring initialization" <<endl;
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    //cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      //cout<<"received radar"<<endl;


      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      /*
      Radar data gives polar coordinates based system, it's not sufficient to determin 
      velocity (vx, vy). but rho (range), phi (bearing) can be used to determin position (px, py)
      Note: Cartesian coordindates for laser: x direction is in direction of motion
       ; y is to-the-left direction
      */
      double m_rho = measurement_pack.raw_measurements_[0];
      double m_phi = measurement_pack.raw_measurements_[1];
      
		  ekf_.x_ << m_rho * cos(m_phi), - m_rho * sin(m_phi), 0, 0;

		  previous_timestamp_ = measurement_pack.timestamp_;

      //cout<< "x_ initialized per laser input:"<< ekf_.x_<<endl;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      //cout<<"received laser"<<endl;
      /**
      Initialize state.
      */
     	//set the state with the initial location and zero velocity
		  ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

		  previous_timestamp_ = measurement_pack.timestamp_;

      //cout<< "x_ initialized per laser input:"<< ekf_.x_<<endl;
  
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;

    //cout << "leaving initialization" <<endl;

    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //cout << "compute elapsed time: " << endl;
  //compute the time elapsed between the current and previous measurements
  //dt - expressed in seconds
	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	

  // TODO (ME):  If two sensor measurements come in simultaneously, the time step
  // between the first measurement and the second measurement would be (close to )
  // zero.
  // process the first arrived as usual; then predict another line again : 
  // predict and update

	previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

	//Modify the F matrix so that the time is integrated
	ekf_.F_(0, 2) = dt;
	ekf_.F_(1, 3) = dt;

	//set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
	            0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
	            dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
	            0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  //cout << ekf_.Q_ << endl;

  //cout<<"before predict()"<<endl;
  ekf_.Predict();
  //cout<<"after predict()" <<endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //cout<<"before radar update"<<endl;

    ekf_.R_ = R_radar_;

    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    //cout<<"after radar update"<<endl;

  } else {
    //cout<<"before laser update"<<endl;

    ekf_.R_ = R_laser_;

    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
    
    //cout<<"after laser update"<<endl;

  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
