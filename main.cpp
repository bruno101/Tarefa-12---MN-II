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

double somaDosQuadradosDosTermosAbaixoDaDiagonal(vector<vector<double>> A, double n) {

  double s = 0;
  for (int i = 0; i < n; i++) {

    for (int j = 0; j < i; j++) {
      s += A[i][j]*A[i][j];
    }

  }
  return s;

}

vector<vector<double>> transpose(vector<vector<double>> A, double n) {

  vector<vector<double>> C(n, vector<double>(n));

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      C[i][j] = A[j][i];
    }
  }

  return C;

}

vector<vector<double>> metodoH_i(vector<vector<double>> A, double i, double n) {

  vector<vector<double>> H(n, vector<double>(n));

  vector<double> w(n), w1(n), N(n), n1(n), e(n);

  double lw = 0;
  for (int k = 0; k < n; k++) {
    H[k][k] = 1;
  }
  for (int k = i+1; k < n; k++) {
    w[k] = A[k][i];
    lw += w[k]*w[k];
  }
  lw = pow(lw,0.5);
  w1[i+1] = lw;

  double s = 0;
  for (int k = 0; k < n; k++) {
    N[k] = w[k] - w1[k];
    s += N[k]*N[k];
  }
  s = pow(s, 0.5);

  for (int k = 0; k < n; k++) {
    n1[k] = N[k]/s;
  }

  for (int j = 0; j < n; j++) {
    for (int k = 0; k < n; k++) {
      H[j][k] = H[j][k] - 2*n1[j]*n1[k];
    }
  }

  return H;

}

tuple<vector<vector<double>>, vector<vector<double>>> metodoH(double n, vector<vector<double>> A) {

  vector<vector<double>> H(n, vector<double>(n)), Hi(n, vector<double>(n)), R(n, vector<double>(n)), A_nova(n, vector<double>(n)), A_velha(n, vector<double>(n)), A_r(n, vector<double>(n));

  for (int i = 0; i < n; i++) {
    H[i][i] = 1;
  }
  A_velha = A;

  for (int i = 0; i < n-2; i++) {

    Hi = metodoH_i(A_velha, i, n);
    A_nova = mult(A_velha, Hi, n);
    A_nova = mult(transpose(Hi, n), A_nova, n);
    A_velha = A_nova;
    H = mult(H, Hi, n);

  }

  return make_tuple(A_nova, H);

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

  for (int j = 0; j < n-1; j++) {

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

void metodoQR(double n, vector<vector<double>> A, double epsilon, double t) {

  vector<vector<double>> P(n, vector<double>(n)), Q(n, vector<double>(n)), R(n, vector<double>(n)), A_nova(n, vector<double>(n)), A_velha(n, vector<double>(n)), A_r(n, vector<double>(n)), H(n, vector<double>(n));
  double val = 100;


  if (t) {
    auto res = metodoH(n, A);
    A_velha = get<0>(res);
    H = get<1>(res);

    cout << "Matriz A calculada pelo metodo de Householder:\n";
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cout << A_velha[i][j] << " | ";
      }
      cout << "\b\b  \n";
    }
    cout << "\n";

  } else {
    A_velha = A;
  }

  for (int i = 0; i < n; i++) {
    P[i][i] = 1;
  }

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
        cout << A_nova[i][j] << " | ";
      }
      cout << "\b\b  \n";
    }
    cout << "\n";

  }

  if (t) {

    A_r = A_nova;
    cout << "Autovetores encontrados a partir das colunas de P:\n";
    for (int i = 0; i < n; i++) {

      cout << "(";
      for (int k = 0; k < n; k++) {
        cout << P[k][i] << ", ";
      }
      cout << "\b\b)\n";
    }
    cout << "\n";

    P = mult(H,P,n);

    A_r = A_nova;
    cout << "Autovalores, e autovetores encontrados a partir das colunas de HP:\n";
    for (int i = 0; i < n; i++) {

      cout << A_r[i][i] << ", (";
      for (int k = 0; k < n; k++) {
        cout << P[k][i] << ", ";
      }
      cout << "\b\b)\n";
    }
    cout << "\n";

  } else {

    A_r = A_nova;
    cout << "Autovalores e autovetores encontrados:\n";
    for (int i = 0; i < n; i++) {

      cout << A_r[i][i] << ", (";
      for (int k = 0; k < n; k++) {
        cout << P[k][i] << ", ";
      }
      cout << "\b\b)\n";
    }
    cout << "\n";

  }

}

int main() {

  int n;

  cout << "Digite a dimensão da matriz: ";
  cin >> n;

  vector<vector<double>> A(n, vector<double>(n)), H(n, vector<double>(n));
  double epsilon;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cout << "Digite o valor de A[" << i+1 << "][" << j+1 << "]: ";
      cin >> A[i][j];
    }
  }

  cout << "Digite o valor da precisão desejada: ";
  cin >> epsilon;

  int t;
  cout << "Usar método de Householder?\n0 - Não\n1 - Sim\n";
  cin >> t;
  cout << "\n";

  metodoQR(n, A, epsilon, t);

}