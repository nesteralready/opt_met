//
// Created by Ilya Nesterenko on 17.04.2021.
//
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

double NormaV(double *V, int n) {
  double norma;
  double s = 0;
  for (int i = 0; i < n; i++)
    s += V[i] * V[i];
  norma = sqrt(s);
  return norma;
}
double Scal(double *v1, double *v2, int n) {
  double s = 0.0;
  for (int i = 0; i < n; i++) {
    s = v1[i] * v2[i] + s;
  }
  return s;
}
double *MatrixByVector(double **H, double *V, int n) {
  double *tmp = new double[n];
  double sum = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      sum += V[j] * H[i][j];
    tmp[i] = sum;
    sum = 0;
  }
  return tmp;
}



double *SoprGrad(double **Q, double * b, int n) {
  int k = 0;
  double *y = new double[n + 1], *x = new double[n + 1],//
  *ax = new double[n + 1], *Qw = new double[n + 1],//ax=Qx
  *z = new double[n + 1], *z1 = new double[n + 1], //z=Qx-b
  *p = new double[n + 1];
  double a, c, nz;

  //X0==
  cout << endl << "Vvedite nachal'noe priblizhenie dlya matrici: " << endl;
  for (int i = 0; i < n; i++)
    cin >> x[i];

  double e;
  cout << "Vvedite to4nost':";
  cin >> e;

  ax = MatrixByVector(Q, x, n);
  for (int i = 0; i < n; i++) z[i] = b[i] - ax[i];

  if (Scal(z, z, n) != 0.0) {
    for (int i = 0; i < n; i++) p[i] = z[i];
    nz = 1000.;
    while (nz > e) {
      Qw = MatrixByVector(Q, p, n);
      a = Scal(z, p, n) / Scal(z, Qw, n);
      for (int i = 0; i < n; i++) {
        y[i] = x[i] + a * p[i];
        z1[i] = z[i] - a * Qw[i];
      }
      nz = NormaV(z1, n);
      c = Scal(z1, Qw, n) / Scal(p, Qw, n);
      for (int i = 0; i < n; i++) {
        p[i] = z1[i] - c * p[i];
        z[i] = z1[i];
        x[i] = y[i];
      }
      k++;
    }
    cout << "k=" << k << endl;
    return y;
  } else {
    cout << "k==" << 1 << endl;
    return x;
  }
}

