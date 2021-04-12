//
// Created by Ilya Nesterenko on 12.04.2021.
//
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>
#include <cassert>
#include <algorithm>


class QuadraticTask {
public:

  QuadraticTask() = default;
  QuadraticTask(unsigned int n_, double k_);

  double eps = 1e-5;
  double eps_lin_search = 1e-5;
  void run(std::vector<double> x_0, Lin_search search){

    fillForm(quadraticForm);

    std::cout << "k = " << k << "   n = " << n << std::endl;
    getForm();
    std::vector res = grad_method(std::move(x_0),search);
    for(auto r : res){
      std::cout << r << " ";
    }
  }

private:
  unsigned int n; //  dimension
  double k; // number obuslovlennosti
  std::vector<double> quadraticForm;

  //sobstvenii  4isla hessiana funkcii
  double lMax;
  double lMin;

  void fillForm(std::vector<double> &form) {
    lMin = 1;
    lMax = k * lMin;
    for (int i = 0; i < n; i++) {
      if (i == 0) {
        form.push_back(lMax);
      } else if (i == n - 1) {
        form.push_back(lMin);
      } else {
        form.push_back(1 + rand() % (int) lMax);
      }
    }
  }

  void getForm() const {
    std::cout << "Current Random Form = " << std::endl;
    for (auto row : quadraticForm) {
        std::cout << row << " ";
      }
    std::cout << std::endl;
  }

  std::vector<double> gradient(std::vector<double> form, std::vector<double> x){

    std::vector<double> gradCoeffs(n);
    for(int i = 0; i < form.size(); i++){
      gradCoeffs[i] = 2 * form[i] * x[i];
    }
    return  gradCoeffs;
  }

  long double f(std::vector<double> x_0) const {
    double res = 0;
    for (int i = 0; i < quadraticForm.size(); i++) {
        res += quadraticForm[i] * pow(x_0[i],2);
    }
    return res;
  }




  [[nodiscard]] std::vector<double> grad_method(std::vector<double> x_0, Lin_search search) {
    int iter = 0;
    std::vector<double> last_point = std::vector<double>(n, INT32_MAX);
    std::vector<double> current_point = std::move(x_0);
    double alpha;
    while(std::abs(f(current_point) - f(last_point)) >= eps){
      last_point = current_point;
      std::vector<double> grad = gradient(quadraticForm, current_point);
      iter++;
    //alpha = golden_sechenie(0,0.5,grad,current_point)[0];
      switch (search) {
      case GOLDEN_SEARCH:
        alpha = golden_sechenie(0,0.5,grad,current_point)[0];
        break;
      case CONST_ALPHA:
        alpha = 0.001;
        break;
      }

      for(int i = 0; i < n; i++) {
        current_point[i] = (current_point[i] - alpha * grad[i]);
      }

    }

    std::cout << iter  << " iter" << std::endl;
    current_point.push_back(f(current_point));
    return current_point;

    // x_k+1 = x_k - a_k* f'(x_k)
    // a_k - shag grad met
    // if f'(x_k) != 0 then  f(x_k+1) < f(x_k)
    //     else x_k - stacion_tochka - minimum
    //
    //vibor shaga
    // 1) naiskoreishego spuska
    //    g_k(a) = f(x_k - a*f'(x_k))
    //   g*_k = inf(g_k(a))
  }

  [[nodiscard]] std::vector<double> golden_sechenie(double a, double b, std::vector<double> grad, std::vector<double> x_k) const {
    assert(a <= b);
    double x_1 = a + (b - a) * (3 - std::sqrt(double(5))) / 2;
    double x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;

    std::vector<double> x1;
    std::vector<double> x2;
    // даже не спрашивайте зачем эти 2 массива, к моменту сдачи  я забуду 99,9%
    std::vector<double> x11(n);
    std::vector<double> x22(n);
    for(int i =0; i < n; i++ ){
      x1.push_back(x_k[i] - x_1 * grad[i]);
      x2.push_back(x_k[i] - x_2 * grad[i]);
    }
    double f_x1 = f(x1);
    double f_x2 = f(x2);
    while ((b - a) > eps_lin_search) {
      if (f_x1 < f_x2) {
        b = x_2;
        x_2 = x_1;
        f_x2 = f_x1;
        x_1 = a + (b - a) * (3 - sqrt(double(5))) / 2;
        for(int i =0; i < n; i++ ){
          x11[i]= (x_k[i] - x_1 * grad[i]);
        }
        f_x1 = f(x11);
      } else {
        a = x_1;
        x_1 = x_2;
        f_x1 = f_x2;
        x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;
        for(int i =0; i < n; i++ ){
          x22[i]= (x_k[i] - x_2 * grad[i]);
        }
        f_x2 = f(x22);
      }
    }
    std::vector<double> end_point;
    for(int i = 0; i < n; i++){
      end_point.push_back((x11[i] + x22[i])/2);
    }
    std::vector<double> res = {(a + b) / 2, f(end_point)};
    return res;
  }



};

QuadraticTask::QuadraticTask(unsigned int n_, double k_) {
  assert(k_ >= 1);
  n = n_;
  k = k_;
}

