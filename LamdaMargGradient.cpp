// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math.hpp>  // pulls in everything from rev/ and prim/
#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::plugins(cpp14)]]
// [[Rcpp::export]]
auto lambdaMargGrad( Eigen::VectorXd theta, int p ) {
  
  // // Declarations //
  // For Stan gradient function:
  double fx;
  Eigen::VectorXd grad_fx;

  // // Compute the gradient using Stan Math //
  stan::math::gradient([ ](Eigen::VectorXd theta, int p) {
    // Assuming theta is organized as a vector (scalar lambda, vector lambda_s, vector Psi_s)
    auto lambda = theta[0];
    // Extract the segment of theta that corresponds to C
    Eigen::VectorXd B = theta.segment(1, p);
    // Extract the segment of theta that corresponds to C
    Eigen::VectorXd C_flat = theta.segment(1 + p, (p * (p + 1)) / 2);
    // Reshape the extracted segment into a p x p matrix
    Eigen::MatrixXd C = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(C_flat.data(), p, p);
    // Define the lambda marginalization function - to be differentiated.
    double result = lambda / stan::math::sqrt(1 + B.transpose() * C * B);
    return result;
  }, theta, fx, grad_fx);
  
  // Return the gradient evaluated at a certain point.
  return grad_fx;}
