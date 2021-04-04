#include <iostream>

#include "OneDimensionSearch.cpp"



int main() {

  OneDimensionSearch lab;
  lab.eps = 1e-5;
  lab.run(0.5,5);

  return 0;

}