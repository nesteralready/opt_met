//
// Created by Ilya Nesterenko on 07.04.2021.
//
#include <cmath>
#include <vector>

using namespace std;
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

  //5 case f(x) = -4x - 2y + 1000x^2 + y^2 + 5
  long double f5(double x, double y) const {
    return 1000*std::pow(x,2) + std::pow(y,2) + 5 - 2 * y - 4 * x;
  }
  long double derivate_f5_x(double x, double y) const {
    return 2000 * x - 4 ;
  }
  long double derivate_f5_y(double x, double y) const {
    return 2 * y - 2;
  }


  //f(x) = sin(1/2x^2 - 1/4y^2 +3)*cos(2x +1+e^y)
  long double f7(double x, double y) const {
    return std::sin(1/2*std::pow(x,2) - 1/4*std::pow(y,2) + 3)* std::cos(2*x + 1 + std::exp(y));
  }
  long double derivate_f7_x(double x, double y) const {
    return (f7(x+0.00001,y) - f7(x,y)) / 0.00001;
  }
  long double derivate_f7_y(double x, double y) const {
    return (f7(x,y+0.00001) - f7(x,y)) / 0.00001;
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


  //f(x) = 100(y-x)^2 + (1-x)^2
// min(1,1) = 0
  long double f2_1(double x, double y) const {
    return 100 * std::pow((y - x),2) + std::pow((1 - x),2);
  }
  long double derivate_f2_1_x(double x, double y) const {
    return 202 * x - 200 * y - 2;
  }
  long double derivate_f2_1_y(double x, double y) const {
    return 200 * (y - x);
  }

  std::vector<std::vector<double>> hessian_f2_1(double x, double y){
    std::vector<std::vector<double>> res{{202,-200},{-200,200}};
    return res;
  }

  //f(x0) = 100(y - x^2)^2 + (1 - x)^2
  //min(1,1) = 0
  long double f3_1(double x, double y) const {
    return 100 * std::pow((y - std::pow(x,2)),2) + std::pow((1 - x),2);
  }
  long double derivate_f3_1_x(double x, double y) const {
    return (f3_1(x + 0.00000001,y) - f3_1(x,y)) / 0.00000001;
  }
  long double derivate_f3_1_y(double x, double y) const {
    return (f3_1(x ,y + 0.00000001) - f3_1(x,y)) / 0.00000001;
  }

  std::vector<std::vector<double>> hessian_f3_1(double x, double y){
    std::vector<std::vector<double>> res{ {1200*std::pow(x,2) - 400 * y + 2,-400 * x},
                                          {-400 * x,200} };
    return res;
  }

  //8 case f(x) = 2^() + 3^()
  long double f4_1(double x, double y) const {
    return  2 * std::exp(-std::pow(((x-1)/2), 2) - std::pow((y-1), 2))
          + 3 * std::exp(-pow(((x-2)/3), 2) - pow((y-3)/2, 2));
  }
  long double derivate_f4_1_x(double x, double y) const {
    return (f4_1(x+0.00001,y) - f4_1(x,y))/0.00001;
  }
  long double derivate_f4_1_y(double x, double y) const {
    return (f4_1(x,y+0.00001) - f4_1(x,y))/0.00001;
  }


  std::vector<std::vector<double>> hessian_f4_1(double x, double y) {
    std::vector<std::vector<double>> res
        {{
             1/2*std::pow((1-x),2) * std::exp(-pow((x - 1), 2) / 4 - pow((y - 1), 2)) -
             2/3*std::pow(3 - y, 2) * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4) -
             std::pow((1-x),2) * std::exp(-pow((x - 1), 2) / 4 - pow((y - 1), 2)) +
             4/27 * std::pow((x-2),2)*std::pow(3 - y, 2) * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4),

             ((4 - 2 * x) / 9) * ((2 * x - 4) / 3) * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4)
              + ((1 - x) / 2) * (x - 1) * std::exp(-pow((x - 1), 2) / 4 - pow((y - 1), 2))
              + 2 * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4) / 3
              + std::exp(-pow((x - 1), 2) / 4 - pow((y - 1), 2)), -400 * x,


          -1/3*(x-2)*(3-y) * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4)
              - 2 * (1 - x) * (y - 1) * std::exp(-0.25 * pow((x - 1), 2) - pow(y - 1, 2))},

          {

          -1/3*(x-2)*(3-y) * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4)
              - 2 * (1 - x) * (y - 1) * std::exp(-0.25 * pow((x - 1), 2) - pow(y - 1, 2)),


          3 / 4 * std::pow(3 - y, 2) * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4)
              - 3 / 2 * std::exp(-pow((x - 2), 2) / 9 - pow((y - 3), 2) / 4) -
              4 * std::exp(-pow((x - 1), 2) / 4 - pow((y - 1), 2))
              + 8 * std::pow((y - 1), 2) * std::exp(-pow((x - 1), 2) / 4 - pow((y - 1), 2))
         }};

    return res;
  }
};

