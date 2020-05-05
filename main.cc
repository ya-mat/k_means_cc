#include <iostream>
#include <Eigen/Dense>

int main()
{
  Eigen::MatrixXd xn;
  Eigen::MatrixXd mu;
  Eigen::VectorXi gk;
  Eigen::VectorXi numk;
  Eigen::MatrixXd::Index minIndex;
  int ni;
  int nn;
  int kk;
  int i;
  int j;
  int l;
  int sum_gk;
  int bef_sum_gk;

  // after :: from file input
  ni = 10;
  kk = 3;

  nn = ni*kk;
  std::cout << "nn " <<  nn << std::endl;

  xn = Eigen::MatrixXd::Random(2, nn);
  gk = Eigen::VectorXi(nn);
  numk = Eigen::VectorXi::Zero(kk);
  std::cout << "xn =" << std::endl << xn << std::endl;

  xn.block(0, nn/kk, 1, nn/kk).array() += 2.0;
  xn.block(0, 2*nn/kk, 1, nn/kk).array() += 4.0;

  mu = xn.block(0, 0, 2, kk);

  std::cout << "mu =" << std::endl << mu << std::endl;

  bef_sum_gk = 0;
  for(i = 0; i < 100; i++){
    std::cout << "loop =" << i << std::endl;
    for(j = 0; j < nn; j++){

// Is which better?
      (mu.colwise() - xn.col(j)).colwise().squaredNorm().minCoeff(&minIndex);
//      (xn.col(j).replicate(1, 3) - mu).colwise().squaredNorm().minCoeff(&minIndex);

      gk(j) = minIndex;
    }
    numk.array() = 0;
    for(l = 0; l < kk; l++){
      mu.col(l).array() = 0.0;
      for(j = 0; j < nn; j++){
	if(gk(j) == l){
	  numk(l) += 1;
	  mu.col(l) += xn.col(j);
	}
      }
// Is which better?
      mu.col(l).array() /= static_cast<double> (numk(l));
//      mu.col(l) /= static_cast<double> (numk(l));
    }
    sum_gk = gk.sum();
    if(bef_sum_gk == sum_gk) break;
    bef_sum_gk = sum_gk;
  }

  std::cout << "numk =" << std::endl << numk << std::endl;
  std::cout << "gk =" << std::endl << gk << std::endl;
  std::cout << "mu =" << std::endl << mu << std::endl;
}
