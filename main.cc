#include <iostream>
#include <Eigen/Dense>
#include <complex>
 
int main()
{
  Eigen::MatrixXd xn;
  Eigen::MatrixXd mu;
  Eigen::VectorXi gk;
  int ni;
  int nn;
  int kk;
  int i;
  int j;
  int l;

  double tmp;
  double now;
  Eigen::VectorXi numk;

  // after :: standard input
  ni = 10;
  kk = 3;

  nn = ni*kk;
  std::cout << "nn " <<  nn << std::endl;

  xn = Eigen::MatrixXd::Random(2, nn);
  gk = Eigen::VectorXi(nn);
  std::cout << "xn =" << std::endl << xn << std::endl;

  for(j = 1; j < kk; j++){
    for(i = 0; i < ni; i++){
      xn(0, j*ni + i) += static_cast<double>(j)*2.0;
    }
  }

  mu = xn.block(0, 0, 2, kk);

  std::cout << "mu =" << std::endl << mu << std::endl;

  for(i = 0; i < 100; i++){
    std::cout << "loop =" << i << std::endl;
    for(j = 0; j < nn; j++){
      now = 0.0;
      for(l = 0; l < kk; l++){
	tmp = std::pow((xn(0, j) - mu(0, l)), 2.0)
	  + std::pow((xn(1, j) - mu(1, l)), 2.0);
	if(now == 0.0){
	  now = tmp;
	  gk(j) = l;
	} else if(tmp <= now){
	  now = tmp;
	  gk(j) = l;
	}
      }
    }
    numk = Eigen::VectorXi::Zero(kk);
    for(l = 0; l < kk; l++){
      mu(0, l) = 0.0;
      mu(1, l) = 0.0;
      for(j = 0; j < nn; j++){
	if(gk(j) == l){
	  numk(l) += 1;
	  mu.col(l) += xn.col(j);
	}
      }
      mu.col(l) /= numk(l);
    }
  }

  std::cout << "xn =" << std::endl << xn << std::endl;
  std::cout << "numk =" << std::endl << numk << std::endl;
  std::cout << "gk =" << std::endl << gk << std::endl;
  std::cout << "mu =" << std::endl << mu << std::endl;


//  xn.block(0, nn/kk, 0, 2*nn/kk-1) = xn.block(0, nn/kk, 0, 2*nn/kk-1) + Eigen::MatrixXd::Constant(0, nn/kk, 1.0);

//  std::cout << xn(1, nn-1) << std::endl;

//  xn.block(0, nn/kk, 0, 2*nn/kk-1) += Eigen::MatrixXd::Constant(1, nn/kk, 1.0);
//  xn.block(0, nn/kk, 0, 2*nn/kk-1) = xn.block(0, nn/kk, 0, 2*nn/kk-1) + Eigen::MatrixXd::Constant(0, nn/kk, 1.0);
//  xn.array.block(0, 2*nn/kk, 0, 3*nn/kk-1) += 2.0

//  xn.col(0)

//  std::cout << "xn =" << std::endl << xn.row(0) << std::endl;
//  std::cout << "xn =" << std::endl << xn.row(1) << std::endl;

//  xn.row(0) = xn.row(0) + Eigen::VectorXd()1.0;
//  xn.row(0) += Eigen::VectorXd::Constant(nn/kk, 1.0);
//  xn.row(0) += Eigen::MatrixXd::Constant(1, nn/kk, 1.0);
//  xn.row(1) = xn.row(1) + 2.0;

//  VectorXcd v(3);
//  v << std::complex<double>(1.0, 2.0), 2, 3;
//  cout << "v =" << endl << v << endl;
//  cout << "m * v =" << endl << m * v << endl;
}

////Boost test
//#include <vector>
//#include <iostream>
//#include <boost/shared_ptr.hpp>
// 
//struct Boo{
//  Boo(){ std::cout<<"Boo "; }
//  ~Boo(){ std::cout<<"~Boo "; }
//};
//typedef boost::shared_ptr<Boo> BooPtr;
// 
//int main(){
//  std::vector<BooPtr> v;
//  v.push_back(BooPtr(new Boo));
//  v.push_back(BooPtr(new Boo));
//}

/*
int main() {
  std::cout << "Hello C++ World" << std::endl;
  return 0;
}
*/
