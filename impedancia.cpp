#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>
#include <random>
#include <iomanip>

using std::complex;

using namespace std;


// protótipo da função da frequência angular
double angular_frequency(double frequency); // protótipo da função

// protótipo da função da impedância real
double impedance_real(int resistor, int capacitor);

// protótipo da função da impedância imaginária
double impedance_image(int resistor, int capacitor);

// random frequency
double generate_data(int f_Hz_initial, int f_Hz_final, float decade, double capacitor, double resistor);


double Gradient(double xo[], double fxo, int n, int MAX_IT, int func, float LAMBDA, float LOWLIMIT, float HIGLIMIT);


double funcao(double S1_c, double S1_r, double S2_c, double S2_r);

// função principal
/****
 * 
 * 
 ****/ 
int main() 
{
  
    int in, final;
    float decade, capacitor, resistor;
    cout << "Digite o valor inicio\n";
    cin >> in;
    
    cout << "Digite o valor final\n";
    cin >> final;
    
    cout << "Digite o valor decade\n";
    cin >> decade;

    cout << "Digite o valor capacitor\n";
    cin >> capacitor;

    cout << "Digite o valor resistor\n";
    cin >> resistor;

    cout << generate_data(in, final, decade, capacitor, resistor);

    return 0;
}

/****
 * 
 * Função angular_frquency()
 * 
 * Descrição: Cálcula a frequência angular.
 * 
 * Parâmetro: frequency (entrada): o valor da frequência.
 * 
 * Retorno: A frequênci angular.
 * 
 ****/
double angular_frequency(double frequency) 
{
    double omega;
    omega = 2*M_PI*frequency;
    
    return omega;
}

/****
 * 
 * Função double impedance_real()
 * 
 * Descrição: Cálcula a impedância real.
 * 
 * Parâmetro: Omega: O valor de entrada para pode cálcular a impedância imaginária.
 * 
 * Retorno: O valor da impedância real.
 * 
 ****/
double impedance_real(double omega)
{
    double z_real, r, c;
    r = 10000;
    c = 1e-8;
    z_real = r/(1+pow(r*c*omega, 2));

    return z_real;   
}

/****
 * 
 * Função double impedance_image()
 * 
 * Descrição: Cálculo a impedância imaginária.
 * 
 * Parâmetro: Omega:  O valor de entrada para pode cálcular a impedância imaginária.
 *
 * Retorno: O valor da impedância imaginária.
 * 
 ****/
double impedance_image(double omega)
{
    double z_image, r, c;   
    r = 10000;
    c = 1e-8;
    z_image = -omega*pow(r, 2)*c/(1 + pow(omega*r*c, 2));
    
    return z_image;
}

/****
 * 
 * Função double generate_data()
 * 
 * Descrição: Gera uma frequência, impedância real e impedância imaginário
 * 
 * Parâmetros: f_Hz_initial(entrada): O valor inicio da frequência. Valor inicial.
 *             f_Hz_final (entrada): O valor final da frequência. Esse valor é o valor final do intervalo.
 *             decade (entrada): Valor por cada década
 * 
 * Retorno: frequência  
 * 
 ****/
double generate_data(int f_Hz_initial, int f_Hz_final, float decade, double capacitor, double resistor)
{
    // 1. Vai gerar um random da frequência e incluir no vetor.
    // 2. Cálcular o omega e colocar em uma vetor também. 
    int n_i, n_f;
    n_i = decade*log10(f_Hz_initial);
    n_f = decade*log10(f_Hz_final);

    
    double f[n_f+1], omega_array[n_f+1], z_real_array[n_f+1], z_image_array[n_f+1];
        
    for (int i = n_i ; i <= n_f+1; i++) {
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(1, -1);

        f[i] = pow(10,  i/decade);
        omega_array[i] = angular_frequency(f[i]);
        z_real_array[i] = (1 + 0.01 * dist(mt)) * impedance_real(omega_array[i]);
        z_image_array[i] = (1 + 0.01 * dist(mt)) * impedance_image(omega_array[i]);
       
    }
   
    ofstream outfile;
    outfile.open("result.txt");
    int i = 0;
    while (i <= n_f) {
        outfile << f[i] << ' ' << z_real_array[i] << ' ' << z_image_array[i] << '\n';
        i++;
    }
    outfile.close();



    return z_real_array[i], z_image_array[i];
}

/****
 * 
 * 
 *
 ****/
double func()
{
    int W;
    double S1_c, S1_r, S2_c, S2_r, r, c, sum_c, sum_r;
    // r -> resitor, c -> capacitor

    S1_c = 2 * W * (r * pow(1 + pow(angular_frequency(i) * r * c, 2), -1) - z_real_array[iters]) * (-2 * pow(r, 3) * pow(angular_frequency(i), 2) * c)/pow(1 + pow(angular_frequency(i) * r * c, 2), 2);
    
    S1_r = 2 * W * (r * pow(1 + pow(angular_frequency(i) * r * c, 2), -1) - z_real_array[iters]) * (-pow(angular_frequency(i) * r * c, 2) + 1)/pow(1 + pow(angular_frequency(i) * r * c, 2), 2);

    S2_c = 2 * W * (-angular_frequency(i) * pow(r,2) * c * pow(1 + angular_frequency(i) + pow(angular_frequency(i) * r * c, 2), -1) - z_image_array(i)) * (-angular_frequency(i) * pow(r, 2) * (pow(-angular_frequency(i) * c * r, 2) + 1))/(pow(1 + pow(angular_frequency(i)* c * r, 2), 2));

    S2_r = 2 * W * (-angular_frequency(i) * pow(r,2) * c * pow(1 + angular_frequency(i) + pow(angular_frequency(i) * r * c, 2), -1) - z_image_array(i)) * (-2 * c * r * angular_frequency(i))/(pow(1 + pow(angular_frequency(i)* c * r, 2), 2));
       

    sum_c = S1_c + S2_c;

    sum_r = S1_r + S2_r;
    
    return sum_c, sum_r;    
}

/****
 * 
 * 
 *
 ****/
#define MAXVAR ??? dimensão da função
double Gradient(double xo[], double fxo, int n, int MAX_IT, int func, float LAMBDA, float LOWLIMIT, float HIGLIMIT)
// x0 é o ponto de partida; fxo é o valor da f.obj em x0; MAX de ITERAÇÕES; parâmetro do método; LIMITES DO DOMÍNIO
{
        double  dir[MAXVAR], xh[MAXVAR], h=0.001F, delta, fx, passo;
        int i, j;

        memcpy(xh,xo, n*sizeof(double));
        for (i=0; i< n; i++){
                xh[i]   = xh[i] + h;
                delta   = funccod[func](xh, n) - fxo;
                dir[i]  = (-1.0F) * delta/h;
                dir[i]  = dir[i] / fabs(dir[i]);
                xh[i]   = xh[i] - h;
        }
        for (i=0; i < MAX_IT; i++){
                for (j=0; j< n; j++){
                        if (dir[j] < 0)
                                passo = LAMBDA*(xh[j] - LOWLIMIT);
                        else
                                passo = LAMBDA*(HIGLIMIT - xh[j]);
                        xh[j] += passo*dir[j];
                }
                fx = funccod[func](xh, n);
         
                if (fx > fxo)
                        break;
                fxo = fx;
        }
        return fxo;
}


