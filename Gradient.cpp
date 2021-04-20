//
// Created by Ilya Nesterenko on 07.04.2021.


struct Point{
  double x;
  double y;
  Point(double x,double y): x(x), y(y){}
  Point() = default;
};

enum Lin_search {
  DIHOTOMIA,
  GOLDEN_SEARCH,
  FIBBONACHI_SEARCH,
  CONST_ALPHA
};

class Gradient{
public:
  double eps = 1e-5;
  double eps_lin_search =  1e-5;
  double delta = eps_lin_search / 2;

  void run(Point x_0, Lin_search search){
    switch (search) {
    case DIHOTOMIA:
      std::cout <<"dihotomia " << std::endl;
      break;
    case GOLDEN_SEARCH:
      std::cout <<"golden se4enie " << std::endl;
      break;
    case FIBBONACHI_SEARCH:
      std::cout <<"fibonacci search " << std::endl;
      break;
    case CONST_ALPHA:
      std::cout <<"const alpha " << std::endl;
      break;
    }

    std::tuple<double,double,double> results = grad_method(x_0,search);
    std::cout << "x= " << std::get<0>(results) << "  y= " << std::get<1>(results)
              << "   z= " << std::get<2>(results) << std::endl ;
    std::cout<< std::endl;
    std::cout<< std::endl;
  }

private:

  Function fun_lib;

  double f(double x, Point grad, Point x_k) const {
    return fun_lib.f4_1(x_k.x - x * grad.x,
                      x_k.y - x * grad.y);
  }


  Point gradient(Point cur){
    double x_grad = fun_lib.derivate_f4_1_x(cur.x,cur.y);
    double y_grad = fun_lib.derivate_f4_1_y(cur.x,cur.y);
    return Point(x_grad,y_grad);
  }


  [[nodiscard]] std::tuple<double, double, double> grad_method(Point x_0, Lin_search search) {
    int iter = 0;
    std::cout << "GRADIENT" << std::endl;
    Point last_point = Point(INT32_MAX,INT32_MAX);
    Point current_point = x_0;
    int gradCounter =0;

    double alpha;
    while(std::abs(fun_lib.f4_1(current_point.x,current_point.y) - fun_lib.f4_1(last_point.x,last_point.y)) >= eps){
      //std::cout << current_point.x << " " << current_point.y <<  std::endl;
      last_point = current_point;
      Point grad = gradient(current_point);
      iter++;
      gradCounter++;
      switch (search) {
      case DIHOTOMIA:
        alpha = std::get<0>(get_min_dihotomia_method(0,0.5,grad,current_point));
        break;
      case GOLDEN_SEARCH:
        alpha = std::get<0>(get_min_golden_sechenie_method(0,0.5,grad,current_point));
        break;
      case FIBBONACHI_SEARCH:
        alpha = std::get<0>(get_min_fibbonaci_method(0,0.5,grad,current_point));
        break;
      case CONST_ALPHA:
        alpha = 0.01;
        break;
      }

      current_point = Point(current_point.x + alpha * grad.x,current_point.y + alpha * grad.y);

     // std::cout << current_point.x << " " << current_point.y << std::endl;

    }

    std::cout << iter  << " iter" << std::endl;
    std::cout << gradCounter  << " grad" << std::endl;
    return std::make_tuple(current_point.x,current_point.y,fun_lib.f4_1(current_point.x,current_point.y));

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


  [[nodiscard]] std::tuple<double, double> get_min_dihotomia_method(double a, double b, Point grad, Point x_k) const {
    assert(a <= b);
    while ((b - a) > eps_lin_search) {
      double x_1 = (a + b) / 2 - delta;
      double x_2 = (a + b) / 2 + delta;
      if (f(x_1,grad,x_k) < f(x_2,grad,x_k)) {
        b = x_1;
      } else {
        a = x_2;
      }
    }
    return std::make_tuple((a + b) / 2, f((a + b) / 2, grad, x_k));
  }

  [[nodiscard]] std::tuple<double, double> get_min_golden_sechenie_method(double a, double b, Point grad, Point x_k) const {
    assert(a <= b);
    double x_1 = a + (b - a) * (3 - std::sqrt(double(5))) / 2;
    double x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;
    double f_x1 = f(x_1, grad, x_k);
    double f_x2 = f(x_2, grad, x_k);

    while ((b - a) > eps_lin_search) {
      if (f_x1 < f_x2) {
        b = x_2;
        x_2 = x_1;
        f_x2 = f_x1;
        x_1 = a + (b - a) * (3 - sqrt(double(5))) / 2;
        f_x1 = f(x_1, grad, x_k);
      } else {
        a = x_1;
        x_1 = x_2;
        f_x1 = f_x2;
        x_2 = a + (b - a) * (sqrt(double(5)) - 1) / 2;
        f_x2 = f(x_2, grad, x_k);
      }
    }
    return std::make_tuple((a + b) / 2, f((a + b) / 2, grad, x_k));
  }


  [[nodiscard]] std::tuple<double, double> get_min_fibbonaci_method(double a, double b, Point grad, Point x_k) const {
    assert(a <= b);

    int realIterationsCounter = 0;
    int functionSolvingCounter = 2;
    int N = 0;
    long double fn1 = 1, fn2 = 1, fn;
    long double target = (b - a) / eps_lin_search;
    while (fn1 < target) {
      fn = fn1 + fn2;
      fn1 = fn2;
      fn2 = fn;
      N++;
    }
    double x1 = a + (double) F(N - 2) / F(N) * (b - a);
    double x2 = a + (double) F(N - 1) / F(N) * (b - a);

    long double xf1 = f(x1, grad, x_k);
    long double xf2 = f(x2, grad, x_k);
    while (realIterationsCounter != N - 2) {
      realIterationsCounter++;
      if (xf1 >= xf2) {
        a = x1;
        x1 = x2;
        xf1 = xf2;
        x2 = a + (double) F(N - realIterationsCounter - 1) / F(N - realIterationsCounter) * (b - a);
        xf2 = f(x2, grad, x_k);

      } else {
        b = x2;
        x2 = x1;
        xf2 = xf1;
        x1 = a + (double) F(N - realIterationsCounter - 2) / F(N - realIterationsCounter) * (b - a);
        xf1 = f(x1, grad, x_k);
      }
      functionSolvingCounter += 1;
    }
    b = a + eps_lin_search;
    return std::make_tuple((a + b) / 2, f((a + b) / 2, grad, x_k));
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