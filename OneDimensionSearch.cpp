//
// Created by Ilya Nesterenko on 04.04.2021.
//

#include <iostream>
#include <cassert>
#include <cmath>
#include <tuple>
#include <chrono>
#include "Function.cpp"
#include <string>

enum Task {
  TASK_ONE,
  TASK_TWO,
  TASK_THREE
};

class OneDimensionSearch {

public:
  double eps = 1e-7;
  double delta = double( eps / 2.0);

  [[nodiscard]] std::tuple<double, double> get_min_dihotomia_method(double a, double b) const {
    assert(a <= b);
    double lenghtOtrezok;
    int functionSolvingCounter = 0;
    int realIterationsCounter = 0;
    double theorIterationCounter = (std::log((b - a) / eps)) / log(2);
    lenghtOtrezok = b - a;
    while ((b - a) > eps) {
      realIterationsCounter = realIterationsCounter + 1;

      std::cout << "i = " << realIterationsCounter << "   a = " << a  << "   b = " << b << std::endl;

      double x_1 = (a + b) / 2 - delta;
      double x_2 = (a + b) / 2 + delta;

      if (f(x_1) < f(x_2)) {
        b = x_1;
        std::cout << "otnoshenie otr = " <<  (b - a) / lenghtOtrezok  << std::endl;
      } else {
        a = x_2;
        std::cout << "otnoshenie otr = " << (b - a) / lenghtOtrezok << std::endl;
      }
      lenghtOtrezok = (b-a);
      functionSolvingCounter += 2;
    }

    // std::cout << "theor number iterations = " << theorIterationCounter << std::endl;
     std::cout << "real number  of iterations = " << realIterationsCounter << std::endl;

     std::cout << "function solving num =  " << functionSolvingCounter << std::endl;

    return std::make_tuple((a + b) / 2, f((a + b) / 2));
  }

  [[nodiscard]] std::tuple<double, double> get_min_golden_sechenie_method(double a, double b) const {

    assert(a <= b);

    int functionSolvingCounter = 2;
    int realIterationsCounter = 0;
    double theorIterationCounter = std::abs(std::log((b - a) / eps) / std::log((std::sqrt(5) - 1) / 2));
    double lenghtOtrezok = b - a;
    double x_1 = a + (b - a) * (3 - std::sqrt(double(5))) / 2;
    double x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;

    double f_x1 = f(x_1);
    double f_x2 = f(x_2);

    while ((b - a) > eps) {

      realIterationsCounter = realIterationsCounter + 1;
       std::cout << "i = " << realIterationsCounter << "   a = " << a  << "   b = " << b << std::endl;
      if (f_x1 < f_x2) {
        b = x_2;
        x_2 = x_1;
        f_x2 = f_x1;
        x_1 = a + (b - a) * (3 - sqrt(double(5))) / 2;

        f_x1 = f(x_1);
        std::cout << "otnoshenie otr = " << (b - a) / lenghtOtrezok   << std::endl;
      } else {
        a = x_1;
        x_1 = x_2;
        f_x1 = f_x2;
        x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;
        f_x2 = f(x_2);
        std::cout << "otnoshenie otr = " << (b - a) / lenghtOtrezok << std::endl;
      }
      lenghtOtrezok = (a-b);
      functionSolvingCounter += 1;
    }
    //std::cout <<  "theor number iterations = " << theorIterationCounter << std::endl;
    std::cout << "real number  of iterations = " << realIterationsCounter << std::endl;

    std::cout << "function solving num =  " << functionSolvingCounter << std::endl;
    return std::make_tuple((a + b) / 2, f((a + b) / 2));
  }


  [[nodiscard]] std::tuple<double, double> get_min_fibbonaci_method(double a, double b) const {
    assert(a <= b);

    int realIterationsCounter = 0;
    int functionSolvingCounter = 2;
    int N = 0;
    long double fn1 = 1, fn2 = 1, fn;
    long double target = (b - a) / eps;
    while (fn1 < target) {
      fn = fn1 + fn2;
      fn1 = fn2;
      fn2 = fn;
      N++;
    }
    //0.61
    //0.38
    double x1 = a + (double) F(N - 2) / F(N) * (b - a);
    double x2 = a + (double) F(N - 1) / F(N) * (b - a);
     double lenghtOtrezok = b - a;
    long double xf1 = f(x1);
    long double xf2 = f(x2);
    while (realIterationsCounter != N - 2) {

      std::cout << "i = " << realIterationsCounter << "   a = " << a  << "   b = " << b << std::endl;
      realIterationsCounter++;
      if (xf1 >= xf2) {
        a = x1;
        x1 = x2;
        xf1 = xf2;
        x2 = a + (double) F(N - realIterationsCounter - 1) / F(N - realIterationsCounter) * (b - a);
        xf2 = f(x2);
        std::cout << "otnoshenie otr = " << (b - a) / lenghtOtrezok << std::endl;

      } else {
        b = x2;
        x2 = x1;
        xf2 = xf1;
        x1 = a + (double) F(N - realIterationsCounter - 2) / F(N - realIterationsCounter) * (b - a);
        xf1 = f(x1);
        std::cout << "otnoshenie otr = " << (b - a) / lenghtOtrezok << std::endl;
      }
      functionSolvingCounter += 1;
    }
    b = a + eps;
    std::cout << "number of  iterations = " << realIterationsCounter << std::endl;

    std::cout << "function solving num =  " << functionSolvingCounter << std::endl;
    return std::make_tuple((a + b) / 2, f((a + b) / 2));
  }


  void run(double a, double b){

    std::cout << "Dihotomia" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::tuple<double, double> resultsDihotomia = get_min_dihotomia_method(a, b);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
    std::cout << "x = " << std::get<0>(resultsDihotomia) << std::endl;
    std::cout << "y = " << std::get<1>(resultsDihotomia) << std::endl;
    std::cout << duration << std::endl;
    std::cout << std::endl;

    std::cout << "GoldenSe4enie" << std::endl;
    auto start_time2 = std::chrono::high_resolution_clock::now();
    std::tuple<double, double> resultsGoldenSe4enie = get_min_golden_sechenie_method(a, b);
    auto end_time2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time2 - start_time2).count();
    std::cout << "x = " << std::get<0>(resultsGoldenSe4enie) << std::endl;
    std::cout << "y = " << std::get<1>(resultsGoldenSe4enie) << std::endl;
    std::cout << duration2 << std::endl;
    std::cout << std::endl;

    std::cout << "Fibbonaci" << std::endl;
    auto start_time3 = std::chrono::high_resolution_clock::now();
    std::tuple<double, double> fibbonaciResults = get_min_fibbonaci_method(a, b);
    auto end_time3 = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time3 - start_time3).count();
    std::cout << "x = " << std::get<0>(fibbonaciResults) << std::endl;
    std::cout << "y = " << std::get<1>(fibbonaciResults) << std::endl;
    std::cout << duration3 << std::endl;

  }

private:

  Function functionLib;

  long double f(double x) const {
      return functionLib.f(x);
  }


  //Fibonacci
  constexpr unsigned long long F(int n) const {
    unsigned long long f1 = 1, f2 = 1, m = 0;
    while (m < n) {
      unsigned long long f = f1 + f2;
      f1 = f2;
      f2 = f;
      m++;
    }
    return f1;
  }

};

