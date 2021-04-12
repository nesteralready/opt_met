#include <iostream>

#include "OneDimensionSearch.cpp"
#include "Gradient.cpp"
#include "QuadraticTask.cpp"

int main() {

  //1 point tests;
  OneDimensionSearch lab;
  lab.eps = 1e-2;
  lab.run(-5,5);



  //2 point tests;
  Gradient grad_test;
  grad_test.eps_lin_search = 1e-5;
  grad_test.eps = 1e-5;
  Point x_0 = Point( 10,10);

  grad_test.run(x_0,DIHOTOMIA);
  grad_test.run(x_0,GOLDEN_SEARCH);
  grad_test.run(x_0,FIBBONACHI_SEARCH);
  grad_test.run(x_0,CONST_ALPHA);

//  //4 point
//  int n = 100;
//  QuadraticTask task(n,10);
//  std::vector<double> x_0;
//  for(int i =0; i< n; i++) {
//    x_0.push_back(1);
//  }
// //= {1,12,1,221,1,1,122,-5,1,12,1,221,1,1,122,-5,1,12,1,221,1,1,122,-5};
//  task.eps = 1e-10;
// // task.run(x_0,CONST_ALPHA);
// task.run(x_0,GOLDEN_SEARCH);
  return 0;

}