 
// teste

// gcc -L/usr/local/lib example.o -lgsl -lgslcblas -lm
//  /usr/lib/x86_64-linux-gnu/libgsl.a
//  /usr/lib/x86_64-linux-gnu/libgsl.so
//  /usr/include/gsl/gsl_sf_legendre.h

//
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

//#include "/usr/include/gsl/gsl_sf_legendre.h"


using namespace std;

/*
double fW(double th, double ph, vector<vector<complex <double>>> Wkq)
{
  double PI = 3.1415926;
  if(th>180.){th=360.-th; ph=ph+180.;}
  double phi=PI*ph/180.; // phi em radianos
  double x=cos(PI*th/180.);
  complex <double> Wc(0.,0.);
  for(int k=0;k<=4;k+=2)
  {
    for(int q=-k;q<=k;q++)
    {
      int m=q;
      double re,im;
      if(m<0){re=cos(m*(PI-phi)); im=-sin(m*(PI-phi));}
      else {re=cos(m*phi); im=-sin(m*phi);}
      complex <double> cf(re,im);
      Wc += Wkq[k][q]*gsl_sf_legendre_sphPlm(k,abs(m),x)*cf;
    }
  }
  return Wc.real();   
}
*/

// Classe das amplitudes dos harmônicos esféricos
class WKQ
{
public:
  WKQ(int km);
  ~WKQ();

  void            set(int, int);
  complex<double> get(int,int);
  double          getReal(int,int);
  double          getIm(int,int);
  
private:
  int kmax;
  complex<double> **Wkq;
}

WKD::WKQ(int km)
{
  kmax = km+1;
  *Wkq = new *complex<double>[kmax];    // k:  de 0 a k
  Wkq  = new complex<double>[2*kmax+1]; // q: -k <= q <= +k
  for(int i = 0; i<kMax; i++)
    for(int j = 0; j<(2*i+1); j++)
      Wkq[i][j] = complex(0.0,0.0);
}

void WKQ::set(int k, int q, double re, double im)
{
  if(k>=0 && k<=kmax && q>=-k && q<=k)
    Wkq[k][kmax+q] = complex(re,im);
  else cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
  return;
}

complex<double> WKQ::get(int k, int q)
{
  if(k>=0 && k<=kmax && q>=-k && q<=k)
    return Wkq[k][k+q];
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

double WKQ::getReal(int k, int q)
{
  if(k>=0 && k<=kmax && q>=-k && q<=k)
    return Wkq[k][k+q].real();
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

double WKQ::getIm(int k, int q)
{
  if(k>=0 && k<=kmax && q>=-k && q<=k)
    return Wkq[k][k+q].im();
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}



int main()
{
  
  int kmax = 4;
  cout << "Define a matriz Wkq" << endl;
  WKQ Wkq(kmax);
  
  cout << " Imprime a matriz" << endl;
  for(int k = 0; k<kmax; k++)
    for(int q = -k; q<=k; q++)
      cout << " Wkq(" << k << "," << q <<") = " <<  Wkq.get(k,q) << endl;
  
//  Wkq.resize(kmax, vector<complex<double>>(kmax));  // zera os elementos
  
//  vector<complex <double>> **Wkq = new *vector<complex <double>>[kmax];
//  Wkq = new vector<complex <double>>[kmax];

//  for(int i = 0; i<kmax; i++)
//    for(int j = 0; j<kmax; j++)
//      Wkq(i,j) = (0.0,0.0);

    cout << "Define alguns valores da matriz Wkq" << endl;
  // valores do Zero -----------------
  Wkq[0][0] = complex<double>(1.406664e-02,0.000000e+00);
//  Wkq[0][0].real() = 1.406664e-02; Wkq[0][0].im() = 0.000000e+00;

  Wkq[2][-2] = complex<double>(-1.026883e-03,-1.476005e-03);
  Wkq[2][-1] = complex<double>(-8.459553e-04,1.611561e-03);
  Wkq[2][0] = complex<double>(-7.571732e-04,-0.000000e+00);
  Wkq[2][1] = complex<double>(8.459553e-04,1.611561e-03);
  Wkq[2][2] = complex<double>(-1.026883e-03,1.476005e-03);

  Wkq[4][-4] = complex<double>(5.458415e-04,-1.472717e-03);
  Wkq[4][-3] = complex<double>(-2.222219e-03,2.716418e-04);
  Wkq[4][-2] = complex<double>(3.012081e-04,4.644176e-04);
  Wkq[4][-1] = complex<double>(-7.038410e-04,1.356695e-03);
  Wkq[4][0] = complex<double>(7.413300e-04,0.000000e+00);
  Wkq[4][1] = complex<double>(7.038410e-04,1.356695e-03);
  Wkq[4][2] = complex<double>(3.012081e-04,-4.644176e-04);
  Wkq[4][3] = complex<double>(2.222219e-03,2.716418e-04);
  Wkq[4][4] = complex<double>(5.458415e-04,1.472717e-03);
  // ---------------------------------

  cout << " Imprime a matriz" << endl;
  for(int i = 0; i<kmax; i++)
    for(int j = 0; j<=i; j++)
      cout << " Wkq(" << i << "," << j<<") = " <<  Wkq[i][j] << endl;
  
  
  
  
  double theta,phi;
  cout << " Gera harmonicos esfericos" << endl;
  
  theta = phi = 0;
  while( theta>=0.0 && phi >=0.0 )
  {
    cout << " Entre com theta:  ";
    cin >> theta;
    cout << " Entre com phi:    ";
    cin >> phi;
  
//    cout << " Wc: " << fW(theta,phi,Wkq) << endl;
  }
  
  
/*
  string nome;  
  cout << "Entre com um nome para o arquivo de saída:";
  cin >> nome;
  ofstream saida;
  saida.open(nome,ios::out);
  
  saida << texto << endl;
  
  saida.close();
*/
  return 0;
}

