#include <iostream>

#include "OneDimensionSearch.cpp"
#include "Gradient.cpp"
#include "QuadraticTask.cpp"
#include "lab2/NewtonMethod.cpp"
#include "lab2/MatrixAPI.cpp"
#include "lab2/SoprGrad.cpp"

#include "lab3/Simplex.cpp"
int main() {
//
//  //1 point tests;
//  OneDimensionSearch lab;
////  lab.eps = 1e-1;
////  lab.run(-5,5);
//  NewtonMethod lab2;
//  auto a = lab2.newton_method({4,4});
//  auto b = lab2.conjGrad({7,9});
//  std::cout <<std::endl;
//  Gradient lab3;
//  lab3.eps = 1e-5;
//


  //lab3.run({-1,3},GOLDEN_SEARCH);


//  //below is about conj grad method
//  int n = 2;
//  auto **Q = new double *[n];
//
//  Q[0] = new double[n]{101, -200};
//  Q[1] = new double[n]{0, 100};
//
//  auto *b = new double[n]{-2, 0};
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




 //lab3 simplex method

// table - last row f(x), first collumn - numbers without variables,

      //for example f(x) = ax1+bx2 -> max, then in last row {0,-a,-b} if ->min then {0,a,b}
      std::vector<std::vector<double>> table = { {7,1,2},
                          {8,2,1},
                          {3,0,1},
                          { 0, -3, -2} };

      std::vector<double> result(2);
      std::vector<std::vector<double>> table_result;

      Simplex S = Simplex(table);

      table_result = S.Calculate(result);
      std::cout << "Решенная симплекс-таблица:" << std::endl;
        for(auto i : table_result){
          for(auto j : i){
            std::cout << j << " ";
          }
          std::cout << std::endl;
        }
       std::cout << std::endl;
       std::cout << "Решение:" << std::endl;
       std::cout << "X[1] = " << result[0] << "   X[2] = "  << result[1];
      std::cout << std::endl;


  return 0;

}