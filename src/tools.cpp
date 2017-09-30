#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() == 0){
    cout << "Error: no estimations were given for rmse!" << endl;
    return rmse;
  }

  if(estimations.size() != ground_truth.size()){
    cout << "Error: estimation and ground truth are different sizes!" << endl;
    return rmse;
  }

  for(int i=0; i < estimations.size(); i++){
    VectorXd diff(estimations.size());
    diff = estimations[i] - ground_truth[i];
    VectorXd residuals = diff.array()*diff.array();
    rmse += residuals;
  }

  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}