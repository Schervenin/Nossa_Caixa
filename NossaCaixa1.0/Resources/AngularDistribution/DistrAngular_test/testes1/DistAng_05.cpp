 
// teste

// g++ -lgsl DistAng_03.cpp -o DistAngNorm20degSolA_03


//  /usr/lib/x86_64-linux-gnu/libgsl.a
//  /usr/lib/x86_64-linux-gnu/libgsl.so
//  /usr/include/gsl/gsl_sf_legendre.h

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;


// Classe das amplitudes dos harmônicos esféricos (Zero)
class WKQ
{
public:
  WKQ(int km);
  ~WKQ();

  void            set(int, int, double, double);
  complex<double> get(int,int);
  double          getReal(int,int);
  double          getIm(int,int);
  
private:
  int kmax;
  complex<double> **Wkq;
};

WKQ::WKQ(int km)
{
  kmax = km+1;
  Wkq = new complex<double>*[kmax];    // k:  de 0 a k
  for(int i=0; i<2*kmax+1; i++)
    Wkq[i]  = new complex<double>[i+1]; // q: -k <= q <= +k
    
  for(int i = 0; i<kmax; i++)
    for(int j = 0; j<(2*i+1); j++)
      Wkq[i][j] = complex<double>(0.0,0.0);
}

WKQ::~WKQ(){;}

void WKQ::set(int k, int q, double re, double im)
{
  if(k>=0 && k<=kmax && q>=-k && q<=k)
    Wkq[k][k+q] = complex<double>(re,im);
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
    return Wkq[k][k+q].imag();
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

// Função para gerar os Harmônicos Espéricos
double fW(double th, double ph, WKQ* Wkq)
{
  double PI = 3.1415926;
  if(th>180.){th=360.-th; ph=ph+180.;}
  double phi=PI*ph/180.; // phi em radianos
  double x=cos(PI*th/180.);
  complex<double> Wc(0.,0.);
  for(int k=0;k<=4;k+=2)
  {
    for(int q=-k;q<=k;q++)
    {
      int m=q;
      double re,im;
      if(m<0){re=cos(m*(PI-phi)); im=-sin(m*(PI-phi));}
      else {re=cos(m*phi); im=-sin(m*phi);}
      complex <double> cf(re,im);
      Wc += Wkq->get(k,q)*gsl_sf_legendre_sphPlm(k,abs(m),x)*cf;
    }
  }
  return Wc.real();   
}



int main()
{
  int kmax = 4;
  cout << "Define a matriz Wkq" << endl;
  WKQ Wkq(kmax);
  
  cout << " Imprime a matriz" << endl;
  for(int k = 0; k<=kmax; k++)
    for(int q = -k; q<=k; q++)
      cout << " Wkq(" << k << "," << q <<") = " <<  Wkq.get(k,q) << endl;
  

  cout << "Define alguns valores da matriz Wkq" << endl;
  // valores do Zero -----------------
  Wkq.set(0,0,1.406664e-02,0.000000e+00);

  Wkq.set(2,-2,-1.026883e-03,-1.476005e-03);
  Wkq.set(2,-1,-8.459553e-04, 1.611561e-03);
  Wkq.set(2, 0,-7.571732e-04,-0.000000e+00);
  Wkq.set(2, 1, 8.459553e-04, 1.611561e-03);
  Wkq.set(2, 2,-1.026883e-03, 1.476005e-03);

  Wkq.set(4,-4, 5.458415e-04,-1.472717e-03);
  Wkq.set(4,-3,-2.222219e-03, 2.716418e-04);
  Wkq.set(4,-2, 3.012081e-04, 4.644176e-04);
  Wkq.set(4,-1,-7.038410e-04, 1.356695e-03);
  Wkq.set(4, 0, 7.413300e-04, 0.000000e+00);
  Wkq.set(4, 1, 7.038410e-04, 1.356695e-03);
  Wkq.set(4, 2, 3.012081e-04,-4.644176e-04);
  Wkq.set(4, 3, 2.222219e-03, 2.716418e-04);
  Wkq.set(4, 4, 5.458415e-04, 1.472717e-03);
  // ---------------------------------

  cout << " Imprime a matriz" << endl;
  for(int k = 0; k<=kmax; k++)
    for(int q = -k; q<=k; q++)
      cout << " Wkq(" << k << "," << q <<") = " <<  Wkq.get(k,q) << endl;
  
  

  double theta,phi;
  cout << " Gera distribuição angular" << endl;
  
  theta = phi = 0;

  while( theta>=0.0 && phi >=0.0 )
  {
    cout << " Entre com theta:  ";
    cin >> theta;
    cout << " Entre com phi:    ";
    cin >> phi;
    cout << " Wc:      " << fW(theta,phi,&Wkq) << endl;
  }
  
  // Grava matriz em Theta - Phi e histogramas das projeções
  // em Theta e Phi
  string nome, nomeTP, nomeTheta, nomePhi;  
  cout << "Entre com um nome para o arquivo de saída:";
  cin >> nome;
  nomeTP    = nome+"TP.dat";
  nomeTheta = nome+"Theta.dat";
  nomePhi   = nome+"Phi.dat";
  ofstream saidaTP, saidaTheta, saidaPhi;
  saidaTP.open(nomeTP,ios::out);
  
  double delta = 2.0;
  
  int nTheta = 180.0/delta+1;
  int nPhi   = 360.0/delta+1;
  
  double* histoTheta = new double[nTheta]; 
  double* histoPhi   = new double[nPhi];
  for(int i = 0; i<nTheta; i++) histoTheta[i] = 0.0;
  for(int i = 0; i<nPhi; i++)   histoPhi[i] = 0.0;
    
  int indiceTheta = 0;
  int indicePhi   = 0;
  
  for(phi = 0.0+delta/2.0; phi < 360.0; phi = phi+delta)
  {
    for(theta = 0.0+delta/2.0; theta < 180.0; theta = theta+delta)
    {
      saidaTP << setw(10) << fW(theta,phi,&Wkq) << " ";
      histoPhi[indicePhi] += fW(theta,phi,&Wkq);  // histograma em phi
    }
    indicePhi++;
    saidaTP << endl;
  }   
  saidaTP.close();
  
  for(theta = 0.0+delta/2.0; theta < 180.0; theta = theta+delta)
  {
    for(phi = 0.0+delta/2.0; phi < 360.0; phi = phi+delta)
      histoTheta[indiceTheta] += fW(theta,phi,&Wkq); // histograma em theta
    indiceTheta++;
  }

  cout << "Grava histograma calculado: " << endl;
  // Grava os histogramas em theta e phi

  saidaPhi.open(nomePhi,ios::out);
  indicePhi = 0;
  for(phi = 0.0+delta/2.0; phi < 360.0; phi = phi+delta)
  {
    saidaPhi << setw(6) << phi << " " << histoPhi[indicePhi] << endl;
    indicePhi++;
  }
  saidaPhi.close();

  
  // Cria histograma de Theta normalizado
  double* histoThetaNorm = new double[nTheta]; 
  for(int i = 0; i<nTheta; i++) histoThetaNorm[i] = 0.0;
  //  calcula área do histoTheta
  double areaHistoTheta = 0.0;
  for(int i = 0; i<nTheta; i++) areaHistoTheta += histoTheta[i];
  for(int i = 0; i<nTheta; i++) histoThetaNorm[i] = histoTheta[i]/areaHistoTheta;

  saidaTheta.open(nomeTheta,ios::out);
  indiceTheta = 0;
  for(theta = 0.0+delta/2.0; theta < 180.0; theta = theta+delta)
  {
    saidaTheta << setw(6) << theta << " " << histoThetaNorm[indiceTheta] << endl;
    indiceTheta++;
  }
  saidaTheta.close();
  
  // Gera valores de Theta por Monte Carlo
  // Inicialização do aleatório
  unsigned long int NT, semente;
  semente = 5168410258;

  gsl_rng_env_setup();

  gsl_rng * ale = gsl_rng_alloc(gsl_rng_mt19937); // Mersenne Twister
  gsl_rng_set(ale, semente);
  
  // Pega o máximo do histograma da distribuição de Theta
  cout << "Procura máximo do histograma de Theta: " << endl;
  double maxDistTheta = -99999999.9;
  for(int i = 0; i<nTheta; i++)
    if(histoThetaNorm[i] > maxDistTheta) maxDistTheta = histoThetaNorm[i];
  cout << "Máximo encontrado: " << maxDistTheta << endl;  

    
  cout << "   Entre com o numero de eventos a gerar: ";
  cin >> NT;
    
  unsigned long int* histoThetaRand = new unsigned long int[nTheta];
  for(int i = 0; i<nTheta; i++) histoThetaRand[i] = 0.0;
  double thetaRand, phiRand;
  for(unsigned long int i = 0; i<NT; i++)
  {
    phiRand   =  360.0*gsl_rng_uniform(ale);  // no intervalo[0,1)
    thetaRand =  180.0*gsl_rng_uniform(ale);  // no intervalo[0,1)
    double yRand = maxDistTheta*gsl_rng_uniform(ale); 
    if(yRand <= fW(thetaRand,phiRand,&Wkq) )  // método da rejeição 
      histoThetaRand[(unsigned long int)round(thetaRand/delta)]++; // incrementa
  }
  
  
  string nomeThetaRand = nome+"ThetaRand.dat";
  ofstream saidaThetaRand;
  saidaThetaRand.open(nomeThetaRand,ios::out);
  indiceTheta = 0;
  for(theta = 0.0+delta/2.0; theta < 180.0; theta = theta+delta)
  {
    saidaThetaRand << setw(6) << theta << " " << histoThetaRand[indiceTheta] << endl;
    indiceTheta++;
  }
  saidaThetaRand.close();
  
  // Calcula número dde rejeitados
  double areaThetaRand = 0.0;
  for(int i = 0; i<nTheta; i++) areaThetaRand += histoThetaRand[i];
  
  unsigned long int Nrej = NT - areaThetaRand;
  cout << " Número de rejeitados: " << Nrej << "  = " << 
       (double)(Nrej)/(double(NT))*100. << "%" << endl;
  


  gsl_rng_free (ale);
  
  cout << "FIM" << endl;
  
  return 0;
}

