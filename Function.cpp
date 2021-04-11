//
// Created by Ilya Nesterenko on 07.04.2021.
//
#include <cmath>
//#include "Gradient.cpp"

class Function {
public:

  Function() = default;


  //1 case x^2 + 2x+ 5 + sin(x)*cos(x)
  long double f(double x) const {
    return std::pow(x, 2) + 2 * x + 5 + std::sin(x) * std::cos(x);
  }
  long double derivate_f(double x) const {
    return 2 * x - std::pow(std::sin(x),2) + std::pow(std::cos(x),2) + 2;
  }

  // 2 case x
  long double f2(double x) const {
    return x;
  }
  long double derivate_f2(double x) const {
    return 1;
  }




  // 3 case f(x) = 10x^2 + y^2
  long double f3(double x, double y) const {
      return 10 * std::pow(x+1,2) + std::pow(y,2) + 5;
  }
  long double derivate_f3_x(double x, double y) const {
      return 20 * x + 20;
  }
  long double derivate_f3_y(double x, double y) const {
    return 2 * y;
  }

  //4 case f(x) = x^2  +2*y^2 + e^(x + y)
  long double f4(double x, double y) const {
    return std::pow(x,2) + 2* std::pow(y,2) + std::exp(x + y);
  }
  long double derivate_f4_x(double x, double y) const {
    return 2 * x +  std::exp(x + y);
  }
  long double derivate_f4_y(double x, double y) const {
    return 4 * y +  std::exp(x + y);
  }

  //5 case f(x) = -4x - 2y + x^2 + y^2 + 5
  long double f5(double x, double y) const {
    return std::pow(x,2) + std::pow(y,2) + 5 - 2 * y - 4 * x;
  }
  long double derivate_f5_x(double x, double y) const {
    return 2 * x - 4 ;
  }
  long double derivate_f5_y(double x, double y) const {
    return 2 * y - 2;
  }

  //6 case f(x) = (1-x)^2 + 100(y-x^2)^2
  long double f6(double x, double y) const {
    return std::pow(1-x,2) + 100* std::pow(y - std::pow(x,2),2);
  }
  long double derivate_f6_x(double x, double y) const {
    return 400 * std::pow(x,3) - 200 * x * y + x -1;
  }
  long double derivate_f6_y(double x, double y) const {
    return 200 * (y - std::pow(x,2));
  }

};
