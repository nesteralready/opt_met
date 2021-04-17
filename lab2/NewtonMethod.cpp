//
// Created by Ilya Nesterenko on 17.04.2021.
//


#include <vector>
#include <iostream>
class NewtonMethod {

public:
  double eps = 1e-10;
  int dim = 2;


  [[nodiscard]] std::vector<double> newton_method(std::vector<double> start_point) {
    std::vector<double> prev_x;
    std::vector<double> cur_x = start_point;
    double alpha;
    do {
      prev_x = cur_x;
      std::vector<double> grad = gradient(cur_x);
      std::vector<std::vector<double>> hess = hessian(cur_x);
      inversion(hess,dim);
      std::vector<double> direction = MatrixByVector(hess,grad,2);

      alpha = (golden_sechenie(0,1,direction,cur_x))[0];
       for(int i = 0; i < dim; i++){
      cur_x[i] = cur_x[i] + alpha * direction[i];
      }
    } while (NormaV({cur_x[0] - prev_x[0], cur_x[1] - prev_x[1]},2 ) >= eps);

    for(int i =0; i < dim; i++){
      std::cout << cur_x[i] << " " << std::endl;
    }
    std::cout << "f(x,y) = " << f(cur_x[0],cur_x[1]) << std::endl;
  }

  void run(std::vector<double> x){

  }


private:

  Function functionLib;

  long double f(double x, double y) const {
    return functionLib.f4_1(x, y);
  }

  [[nodiscard]] std::vector<double> gradient(std::vector<double> x_k){
    std::vector<double> res;
    res.push_back(functionLib.derivate_f4_1_x(x_k[0],x_k[1]));
    res.push_back(functionLib.derivate_f4_1_y(x_k[0],x_k[1]));
    return res;
  }

  [[nodiscard]] std::vector<std::vector<double>> hessian(std::vector<double> x_k){
    std::vector<std::vector<double>> res;
    res = functionLib.hessian_f4_1(x_k[0], x_k[1]);
    return res;
  }

  std::vector<double> MatrixByVector(std::vector<std::vector<double>> H, std::vector<double>  V, int n) {
    std::vector<double> tmp(n);
    double sum = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        sum += V[j] * H[i][j];
      }
      tmp[i] = sum;
      sum = 0;
    }
    return tmp;
  }

  void inversion(std::vector<std::vector<double>> &A, int N)
  {
    long double temp;

   std::vector<std::vector<double>> E;

    for (int i = 0; i < N; i++)
      E.push_back(std::vector<double>(N,0.0));

    for (int i = 0; i < N; i++)
      E[i][i] = 1.0;


    for (int k = 0; k < N; k++)
    {
      temp = A[k][k];
      if( temp  == 0){
        temp += 0.0000001;
      }
      for (int j = 0; j < N; j++)
      {
        A[k][j] /= temp;
        E[k][j] /= temp;
      }

      for (int i = k + 1; i < N; i++)
      {
        temp = A[i][k];
        for (int j = 0; j < N; j++)
        {
          A[i][j] -= A[k][j] * temp;
          E[i][j] -= E[k][j] * temp;
        }
      }
    }

    for (int k = N - 1; k > 0; k--)
    {
      for (int i = k - 1; i >= 0; i--)
      {
        temp = A[i][k];
        for (int j = 0; j < N; j++)
        {
          A[i][j] -= A[k][j] * temp;
          E[i][j] -= E[k][j] * temp;
        }
      }
    }

    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        A[i][j] = E[i][j];

    E.clear();
  }

  [[nodicard]] double NormaV(std::vector<double> V, int n) {
    double norma;
    double s = 0;
    for (int i = 0; i < n; i++)
      s += V[i] * V[i];
    norma = std::sqrt(s);
    return norma;
  }

  long double l_lin_search(double x){

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

  [[nodiscard]] std::vector<double> golden_sechenie(double a, double b, std::vector<double> grad, std::vector<double> x_k) const {
    assert(a <= b);
    double x_1 = a + (b - a) * (3 - std::sqrt(double(5))) / 2;
    double x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;

    std::vector<double> x1;
    std::vector<double> x2;
    // даже не спрашивайте зачем эти 2 массива, к моменту сдачи  я забуду 99,9%
    std::vector<double> x11(dim);
    std::vector<double> x22(dim);
    for(int i =0; i < dim; i++ ){
      x1.push_back(x_k[i] - x_1 * grad[i]);
      x2.push_back(x_k[i] - x_2 * grad[i]);
    }
    double f_x1 = f(x1[0],x1[1]);
    double f_x2 = f(x2[0],x2[1]);
    while ((b - a) > 1e-6) {
      if (f_x1 < f_x2) {
        b = x_2;
        x_2 = x_1;
        f_x2 = f_x1;
        x_1 = a + (b - a) * (3 - sqrt(double(5))) / 2;
        for(int i =0; i < dim; i++ ){
          x11[i]= (x_k[i] - x_1 * grad[i]);
        }
        f_x1 = f(x11[0],x11[1]);
      } else {
        a = x_1;
        x_1 = x_2;
        f_x1 = f_x2;
        x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;
        for(int i =0; i < dim; i++ ){
          x22[i]= (x_k[i] - x_2 * grad[i]);
        }
        f_x2 = f(x22[0],x22[1]);
      }
    }
    std::vector<double> res = {(a + b) / 2};
    return res;
  }
};