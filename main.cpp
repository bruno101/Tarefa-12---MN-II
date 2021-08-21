#include <iostream>
#include <vector>
#include <tuple>
#include <math.h>
using namespace std;

vector<vector<double>> mult(vector<vector<double>> A, vector<vector<double>> B, double n) {

  vector<vector<double>> C(n, vector<double>(n));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = 0;
      for (int k = 0; k < n; k++) {
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }

  return C;
}

vector<double> multV(vector<vector<double>> A, vector<double> v, double n) {

  vector<double> r(n);

  for (int i = 0; i < n; i++) {
    r[i] = 0;
    for (int j = 0; j < n; j++) {
      r[i] += A[i][j]*v[j];
    }
  }

  return r;

}

double somaDosQuadradosDosTermosAbaixoDaDiagonal(vector<vector<double>> A, double n) {
  double s = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j > i; j++) {
      s += A[i][j]*A[i][j];
    }
  }
  return s;
}

vector<vector<double>> matrizJacobiBaseadaNoElemento_ij_DeRvelha(vector<vector<double>> A, int i, int j, int n) {

  vector<vector<double>> Jij(n, vector<double>(n));
  double theta;
  double epsilon = 0.000001;

  for (int i = 0; i < n; i++) {
    Jij[i][i] = 1;
  }

  if (abs(A[i][j]) <= epsilon) {
    return Jij;
  }

  if (abs(A[j][j]) <= epsilon) {
    double PI = 3.14159265358979323846;
    if (A[i][j] < 0) {
      theta = PI/2.0;
    } else {
      theta = -PI/2.0;
    }
  } else {
    theta = atan(-A[i][j]/A[j][j]);
  }

  Jij[i][i] = cos(theta);
  Jij[j][j] = cos(theta);
  Jij[i][j] = sin(theta);
  Jij[j][i] = -sin(theta);

  return Jij;

}

tuple <vector<vector<double>>, vector<vector<double>>> decomposicaoQR(vector<vector<double>> A, double n) {

  vector<vector<double>> QT(n, vector<double>(n)), Jij(n, vector<double>(n)), R_nova(n, vector<double>(n)), R_velha(n, vector<double>(n)), R(n, vector<double>(n)), Q(n, vector<double>(n));

  for (int i = 0; i < n; i++) {
    QT[i][i] = 1;
  }
  R_velha = A;

  for (int j = 0; j < n; j++) {
    for (int i = j+1; i < n; i++) {
      Jij = matrizJacobiBaseadaNoElemento_ij_DeRvelha(R_velha, i, j, n);
      R_nova = mult(Jij, R_velha, n);
      R_velha = R_nova;
      QT = mult(Jij, QT, n);
    }
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Q[j][i] = QT[i][j];
    }
  }
  R = R_nova;

  return make_tuple(Q,R);

} 

void metodoQR(double n, vector<vector<double>> A, vector<double> v0, double epsilon) {

  vector<vector<double>> P(n, vector<double>(n)), Q(n, vector<double>(n)), R(n, vector<double>(n)), A_nova(n, vector<double>(n)), A_velha(n, vector<double>(n)), A_r(n, vector<double>(n));
  vector<double> Lamb;
  double val = 100;

  for (int i = 0; i < n; i++) {
    P[i][i] = 1;
  }
  A_velha = A;

  while (val > epsilon) {

    auto res = decomposicaoQR(A_velha, n);
    Q = get<0>(res); R = get<1>(res);
    A_nova = mult(R,Q,n);
    A_velha = A_nova;
    P = mult(P,Q,n);
    val = somaDosQuadradosDosTermosAbaixoDaDiagonal(A_nova, n);

    cout << "Matriz diagonal calculada nessa iteração:\n";
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cout << A_r[i][j] << " | ";
      }
      cout << "\b\b\n";
    }
    cout << "\n";

  }

  A_r = A_nova;
  cout << "Autovalores e autovetores encontrados:\n";
  for (int i = 0; i < n; i++) {

    cout << A_r[i][i] << ", (";
    vector<double> aV(n);
    for (int k = 0; k < n; k++) {
      aV[k] = A_r[i][i];
    }
    aV = multV(P, aV, n);
    for (int j = 0; j < n; j++) {
      cout << aV[j] << ", ";
    }
    cout << "\b\b)\n";
  }
  cout << "\n";

}

int main() {

  int n;

  cout << "Digite a dimensão da matriz: ";
  cin >> n;

  vector<double> v0(n);
  vector<vector<double>> A(n, vector<double>(n));
  vector<double> v(n);
  double lambda;
  double epsilon;

  for (int i = 0; i < n; i++) {
    v0[i] = 1;
    for (int j = 0; j < n; j++) {
      cout << "Digite o valor de A[" << i+1 << "][" << j+1 << "]: ";
      cin >> A[i][j];
    }
  }

  cout << "Digite o valor da precisão desejada: ";
  cin >> epsilon;

  metodoQR(n, A, v0, epsilon);


}