#include <iostream>

#include "OneDimensionSearch.cpp"
#include "Gradient.cpp"
#include "QuadraticTask.cpp"
#include "lab2/NewtonMethod.cpp"
#include "lab2/MatrixAPI.cpp"
#include "lab2/SoprGrad.cpp"
int main() {
//
//  //1 point tests;
  OneDimensionSearch lab;
//  lab.eps = 1e-1;
//  lab.run(-5,5);
  NewtonMethod lab2;
  auto a = lab2.newton_method({100,100});


  //below is about conj grad method
//  int n = 2;
//  auto **Q = new double *[n];
//
//  Q[0] = new double[n]{3, 4};
//  Q[1] = new double[n]{4, -3};
//
//  auto *b = new double[n]{1, 5};
//
//  double *k = SoprGrad(Q, b, n);
//  for (int i = 0; i < n; i++) {
//    cout << k[i] << "  ";
//  }

//  //2 point tests;
//  Gradient grad_test;
//  grad_test.eps_lin_search = 1e-5;
//  grad_test.eps = 1e-10;
//  Point x_0 = Point( -0.33,0.33);
//
//
//  grad_test.run(x_0,DIHOTOMIA);
//  grad_test.run(x_0,GOLDEN_SEARCH);
//  grad_test.run(x_0,FIBBONACHI_SEARCH);
//  grad_test.run(x_0,CONST_ALPHA);

//  //4 point
//  int n = 100;
//  QuadraticTask task(n,10);
//  std::vector<double> x_0;
//  for(int i =0; i< n; i++) {
//    x_0.push_back(1);
//  }
 //= {1,12,1,221,1,1,122,-5,1,12,1,221,1,1,122,-5,1,12,1,221,1,1,122,-5};
//
//  for(int i =1 ; i < 102;  i += 10){
//    for(int k =1 ; k < 102;  k += 10) {
//      QuadraticTask task(i, k);
//      task.eps = 1e-10;
//      std::vector<double> x_0;
//      for (int j = 0; j < i; j++) {
//        x_0.push_back(1);
//      }
//      task.run(x_0,GOLDEN_SEARCH);
//    }
//  }
 // task.run(x_0,CONST_ALPHA);
  return 0;

}