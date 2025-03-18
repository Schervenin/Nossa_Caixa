 
// teste
// esta versão usa leitura de arquivo

// g++ -lgsl DistAng_xx.cpp -o DistAng_xx


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


// --------- Classe das amplitudes dos harmônicos esféricos (Zero)
class WKQ
{
 public:
   WKQ(int km);
   ~WKQ();

   void            set(int, int, double, double);
   int             getKmax() {return kmax;}
   complex<double> get(int,int);
   double          getReal(int,int);
   double          getIm(int,int);
   void            ReadFromFile(std::string fileName);
   void            ResetValues();
   double          SphHarm(double th, double ph);
  
 private:
   int kmax, indmax;
   complex<double> **Wkq;
};

// ----- Construtor
WKQ::WKQ(int km)
{
  kmax = km;
  indmax = kmax+1;  // indices de 0 a kmax
  Wkq = new complex<double>*[indmax];    // k:  de 0 a k
  for(int i=0; i<indmax; i++)
    Wkq[i]  = new complex<double>[2*(indmax-1)+1]; // q: -k <= q <= +k
    
  ResetValues();  // inicializa com W_00 = (1.0, 0.0), o resto nulo
}

// ----- Destrutor
WKQ::~WKQ(){;}


// ----- Insere valores
void WKQ::ReadFromFile(std::string fileName)
{
  std::ifstream WklFile;
  WklFile.open(fileName,std::ios::in);
  if(!WklFile)
  {
      cout << "Can not open file " << fileName << endl;
      cout << "Only W_00 = (1.0,0.0) will be defined" << endl;
      ResetValues();
      return;
  }
  
  int k,l;
  double realP, imP;
  while( !WklFile.eof() ) 
  {
    WklFile >> k;
    if(k<0)
    {
      cout << " k < 0. Reading of file " << fileName << " will stop." << endl;
      return;
    }
      
    if(k>kmax)
    {
      cout << " k > " << kmax << " in file " << fileName << ". File reading will stop." << endl;
      cout << " Create a new WKQ object with larger kmax if needed. " << endl;
      return;      
    }
      
    cout << "k: " << k << endl;
    for(int i = 0; i < 2*k+1; i++)
    {
      cout << "i: " << i << endl;
        
      if(!(WklFile >> l >> realP >> imP))
      {
        if(WklFile.eof())
        {
          WklFile.close();
          return;
        }
        else
        {
          cout << " WRONG LINES IN FILE " << fileName << endl;
          cout << "   For each k value, 2*k+1 lines must be given" << endl;
          cout << "   Only W_00 = (1.0,0.0) will be defined" << endl;
          ResetValues();
          return;
        }
      }
      set(k,l,realP,imP);
    }    
  } 
  WklFile.close();
}

// ----- Reset
void WKQ::ResetValues()
{
  for(int i = 0; i<indmax; i++)
    for(int j = 0; j<(2*(indmax-1)+1); j++)
      Wkq[i][j] = complex<double>(0.0,0.0);
  Wkq[0][0] = complex<double>(1.0,0.0); // para k=0, m=0
}

// ----- Insere valores
void WKQ::set(int k, int q, double re, double im)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
  {
    Wkq[k][k+q] = complex<double>(re,im);
    cout << "k = " << k << " - q = " << q << " - k+q = " << k+q << endl;
    cout << " Set Wkq("<<k<<","<<q<<")= " <<Wkq[k][k+q]<<endl;
  }
  else cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
  return;
}

// ----- Lê valores
complex<double> WKQ::get(int k, int q)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
    return Wkq[k][k+q];
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

double WKQ::getReal(int k, int q)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
    return Wkq[k][k+q].real();
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

double WKQ::getIm(int k, int q)
{
  if(k>=0 && k<=indmax && q>=-k && q<=k)
    return Wkq[k][k+q].imag();
  else
  {
    cout << "Índices inválidos: k = " << k << "; q = " << q << endl;
    return 0;
  }
}

// Função para gerar os Harmônicos Espéricos
double WKQ::SphHarm(double th, double ph)
{
  double PI = 3.14159265359;
  if(th>PI){th=2.0*PI-th; ph=ph+PI;}
  double x=cos(th);
  complex<double> Wc(0.,0.);
//  for(int k=0;k<=Kmax;k+=2)
  for(int k=0;k<=kmax;k++)    
  {
    for(int q=-k;q<=k;q++)
    {
      double re,im;
      if(q>=0) {re=cos(q*ph); im=-sin(q*ph);}
      else {re=cos(q*(PI-ph)); im=-sin(q*(PI-ph));}
      complex <double> cf(re,im);
      Wc += get(k,q)*gsl_sf_legendre_sphPlm(k,abs(q),x)*cf;
    }
  }
  return Wc.real();   
}
// ---------------------------------------------------------


int main()
{
  int kmax = 4;
  cout << "Define os coeficientes Wkq" << endl;
  WKQ Wkq(kmax);
  
  cout << " Imprime a matriz" << endl;
  for(int k = 0; k<=kmax; k++)
    for(int q = -k; q<=k; q++)
      cout << " Wkq(" << k << "," << q <<") = " <<  Wkq.get(k,q) << endl;
    
  string arquivoWkl;
  cout << " Entre com o nome do arquivo dos Wkl: ";
  cin >> arquivoWkl;
 
  Wkq.ReadFromFile(arquivoWkl);
  
  // ---------------------------------
  cout << " Verificação: imprime a matriz" << endl;
  for(int k = 0; k<=kmax; k++)
    for(int q = -k; q<=k; q++)
      cout << " Wkq(" << k << "," << q <<") = " <<  Wkq.get(k,q) << endl;
  // ---------------------------------
    
  double theta,phi;
  cout << " Gera distribuição angular" << endl;
  
  theta = phi = 0.0;
  double PI    = 3.14159265359;
  double TwoPI = PI*2.0;
  double d2r   = PI/180.0;
  while( theta>=0.0 && phi >=0.0 )
  {
    cout << " Entre com theta (graus):  ";
    cin >> theta;
    cout << " Entre com phi (graus):    ";
    cin >> phi;
    cout << " Wc: " << Wkq.SphHarm(theta*d2r,phi*d2r) << endl;
  }
  
  double delta, deltaGraus;
  cout << " Entre com a variação angular em graus: ";
  cin >> deltaGraus;
  delta = deltaGraus*d2r;
  
  
  // ------------- Grava matriz em Theta - Phi
  string nome, nomeTP, nomeTPMC;
  cout << " Entre com um nome para o arquivo de saída: ";
  cin  >> nome;
  nomeTP = nome+"TP.dat";
  ofstream saidaTP;
  saidaTP.open(nomeTP,ios::out);
  
  int nTheta = ceil(180.0/deltaGraus)+1;
  int nPhi   = ceil(360.0/deltaGraus)+1;
        
  // Guarda os valores máximo e mínimo
  double maxfW = -9.9e10;
  double minfW =  9.9e10;
  double temp  = 0.0;
  for(phi = 0.0+delta/2.0; phi < TwoPI; phi = phi+delta)
  {
    for(theta = 0.0+delta/2.0; theta < PI; theta = theta+delta)
    {
      temp = Wkq.SphHarm(theta,phi);
      saidaTP << setw(10) << temp << " ";
      if(temp > maxfW) maxfW = temp;
      if(temp < minfW) minfW = temp;
    }
    saidaTP << endl;
  }   
  saidaTP.close();
  cout << " mínimo de fW = " << minfW << endl;
  cout << " máximo de fW = " << maxfW << endl;
  
  // --------------------- Gera por Método de Monte Carlo

  // Inicialização do aleatório
  unsigned long int NT, semente;
  semente = 5168410258L;

  gsl_rng_env_setup();

  gsl_rng * ale = gsl_rng_alloc(gsl_rng_mt19937); // Mersenne Twister
  gsl_rng_set(ale, semente);
  
  cout << "   Entre com o numero de eventos a gerar: ";
  cin >> NT;
  
  unsigned long int** DistribRand = new unsigned long int*[nPhi]; // lines
  for(int i=0; i<nPhi; i++)
    DistribRand[i] = new unsigned long int[nTheta]; // columns
  
  for(int i=0; i<nPhi; i++)
    for(int j=0; j<nTheta; j++)
      DistribRand[i][j] = 0;

  int indiceTheta = 0;
  int indicePhi   = 0;
  double phiRand, thetaRand;
    
  for(unsigned long int i = 0; i<NT; i++)
  {
    // método da rejeição 
    phiRand   =  TwoPI*gsl_rng_uniform(ale);  // random no intervalo[0,1)
    thetaRand =  PI*gsl_rng_uniform(ale);  
//    thetaRand =  acos(1.0-2.0*gsl_rng_uniform(ale));  
    double yRand = maxfW*gsl_rng_uniform(ale);
    if(yRand <= minfW)
    {
      int indTh = (int)(round(thetaRand/delta));
      int indPh = (int)(round(phiRand/delta)); 
      DistribRand[indPh][indTh]++;
    }
    else if(yRand <= Wkq.SphHarm(thetaRand,phiRand) )
    {
      int indTh = (int)(round(thetaRand/delta));
      int indPh = (int)(round(phiRand/delta)); 
      DistribRand[indPh][indTh]++;
    }
  }
  
  // Calcula número de rejeitados
  unsigned long int numAccepted = 0L;
  for(int i=0; i<nPhi; i++)
    for(int j=0; j<nTheta; j++)
      numAccepted += DistribRand[i][j];
  
  unsigned long int Nrej = NT - numAccepted;
  cout << " Número de rejeitados: " << Nrej << "  = " << 
       (double)(Nrej)/(double(NT))*100. << "%" << endl;

  
  // Monta matriz da distribuição normalizada
  double** DistribRandNorm = new double*[nPhi]; // lines
  for(int i=0; i<nPhi; i++)
    DistribRandNorm[i] = new double[nTheta]; // columns
  
  for(int i=0; i<nPhi; i++)
    for(int j=0; j<nTheta; j++)
      DistribRandNorm[i][j] = (double)(DistribRand[i][j])/(double)(numAccepted);
  
  // --- Grava em arquivo
  nomeTPMC = nome+"TPMC.dat";
  ofstream saidaTPMC;
  saidaTPMC.open(nomeTPMC,ios::out);
  for(phi = 0.0+delta/2.0; phi < TwoPI; phi = phi+delta)
  {
    for(theta = 0.0+delta/2.0; theta < PI; theta = theta+delta)
    {
      int indTh = (int)(round(theta/delta));
      int indPh = (int)(round(phi/delta)); 
      saidaTPMC << setw(10) << DistribRandNorm[indPh][indTh] << " ";
    }
    saidaTPMC << endl;
  }   
  saidaTPMC.close();
  
  
  gsl_rng_free (ale);

  cout << "FIM" << endl;
  
  return 0;
}

