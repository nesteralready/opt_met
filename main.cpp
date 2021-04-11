#include <iostream>

#include "OneDimensionSearch.cpp"
#include "Gradient.cpp"


int main() {

  //1 point tests;
  OneDimensionSearch lab;
  lab.eps = 1e-2;
  lab.run(-5,5);



  //2 point tests;
  Gradient grad_test;
  grad_test.eps_lin_search = 1e-10;
  grad_test.eps = 1e-7;
  Point x_0 = Point(10,10);

  grad_test.run(x_0,DIHOTOMIA);
  grad_test.run(x_0,GOLDEN_SEARCH);
  grad_test.run(x_0,FIBBONACHI_SEARCH);
  grad_test.run(x_0,CONST_ALPHA);

  return 0;

}