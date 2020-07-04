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



    return 0;
}

/****
 * 
 * 
 *
 ****/
double funcao(double S1_c, double S1_r, double S2_c, double S2_r)
{
    double S1_C,S1_R, capacitor_old,resistor_old, error, c_error, r_error, alpha, S1_C_2, S1_R_2;

    error = 0.001;
    alpha = 0.001;

    int W = 1;

    
    // interações 
    unsigned int iters = 0;
    unsigned int max_iters = 10000;
    c_error = error + 1;
    r_error = error + 1;
    
    while ((error < c_error) && (error < r_error) && (iters < max_iters)) {
        
        capacitor_old = capacitor;
        resistor_old = resistor;
        
        // gradient em função do capacitor        
        S1_C = 2 * W * (resistor * pow(1+pow(omega_array[iters] * resistor * capacitor, 2), -1) - z_real_array[iters]) * (-2 * pow(resistor, 3) * pow(omega_array[iters], 2) * capacitor)/pow(1 + pow(omega_array[iters] * resistor * capacitor, 2), 2);

        // gradient em função do resistor
        // S1_R = 2 * W * (resistor - z_real_array[iters] * (1 + pow(omega_array[iters]*resistor*capacitor, 2))) * (-1 + 3 * pow(omega_array[iters]*resistor*capacitor, 2))/pow(resistor, 2)*pow(1 + pow(omega_array[iters]*resistor*capacitor, 2), 3);
        S1_R = 2 * W * (resistor * pow(1+pow(omega_array[iters] * resistor * capacitor, 2), -1) - z_real_array[iters]) * (-pow(omega_array[iters] * resistor * capacitor, 2) + 1)/pow(1 + pow(omega_array[iters] * resistor * capacitor, 2), 2);

        // capacitor = capacitor - alpha * S1_C
        capacitor -= alpha * S1_C;
        
        resistor -= alpha * S1_R;

        // delta C => 0.001
        // delta R => 0.001
        c_error = abs(capacitor_old - capacitor);
        r_error = abs(resistor_old - resistor);

        // c_error = abs(capacitor - capacitor_old);
        // r_error = abs(resistor - resistor_old);
        iters++;

        cout << "c_erro: " << c_error << "\n";
        cout << "r_error: " << r_error << "\n\n";

        cout << "S1_C: " << S1_C << "\n";
        // cout << "S1_C_2: " << S1_C_2 << "\n\n";

        cout << "S1_R: " << S1_R << "\n";
        // cout << "S1_R_2: " << S1_R_2 << "\n\n";

    return 0;    
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










    
    
    


    }









}

/****
 * 
 * Função ouble derivative_C(double W=1, double alpha=0.01)
 * 
 * Descrição : Função que é responsável pela método de gradiente.
 *             Nessa função tem as derivadas S1_C, S1_R.
 * 
 * Parâmetros: W=1: É o peso da equação que já está definido na entrada
 *             aplha: Taxa de aprendizado
 * 
 * 
 * Retorno: Novo capacitor e o resistor           
 * 
 ****/
// double gradient_descent(double W, double alpha, double resistor, double capacitor, unsigned int max_iters, double omega_array[], double z_real_array[], double z_image_array[])
// {
//     double S1_C,S1_R, capacitor_old,resistor_old, error, c_error, r_error;
//     error = 0.001;

    
//     // interações 
//     unsigned int i = 0;
//     c_error = error + 1;
//     r_error = error + 1;
    
//     while (error < c_error  && error < r_error && i < max_iters) {
        
//         capacitor_old = capacitor;
//         resistor_old = resistor;
        
//         // gradient em função do capacitor
//         S1_C = 2 * W * (2* pow(omega_array[i], 2) * resistor * capacitor * (resistor - z_real_array[i] * (1 + pow(omega_array[i]*resistor*capacitor, 2))))/pow(1 + pow(omega_array[i]*resistor*capacitor, 2), 3);
        
//         // gradient em função do resistor
//         S1_R = 2 * W * (resistor - z_real_array[i] * (1 + pow(omega_array[i]*resistor*capacitor, 2))) * (-1 + 3 * pow(omega_array[i]*resistor*capacitor, 2))/pow(resistor, 2)*pow(1 + pow(omega_array[i]*resistor*capacitor, 2), 3);


//         // capacitor = capacitor - alpha * S1_C
//         capacitor -= alpha * S1_C;
        
//         resistor -= alpha * S1_R;

//         // delta C => 0.001
//         // delta R => 0.001
//         c_error = abs(capacitor_old - capacitor);
//         r_error = abs(resistor_old - resistor);
//         i++;

//         cout << "c_erro: " << c_error << "\n";
//         cout << "r_error: " << r_error << "\n\n";


//     }
//         return 0;
// }



 
 
 
 // gradient_descent(1, 0.001, 0.1, 0.1, 10000);




    // Gradient Z' => dS1/dC
    // int W = 1;
    // double S1_C;
  

    // S1_C => capacitor
    // S1_C = 2 * W * (2* O * resistor * capacitor * (resistor - z_real_array[i] * (1 + O * R * C )))/pow(1+ O * C * R, 3);

    // Gradient Z' => dS1/dR
    // double S1_R;
    // S1_R => resistor
    // S1_R = 2 * W * (resistor - z_real_array[i] * (1 + O * R * C)) * (-1 + 3 * O * R * C)/R*pow(1+O*R*C, 3);

    // Gradient Z" => dS2/dC
    // double S2_C, S2_R;

    // S2_C => capacitor 
    // S2_C = 2 * W * (-omega_array[i] * R * capacitor - z_image_array[i] * (1 + O[i]*R*C))*(1 + 3*O[i]*R*C)/O[i]*R*C*pow(1+O[i]*R*C, 3);


    // S2_R => resistor
    // S2_R = 2 * W * (2 * (-omega_array[i] * R * capacitor - z_image_array[i]*(1 + O[i]*R*C))*(2*O[i]*R*C + 1))/omega_array[i]*pow(resistor, 3)*capacitor*pow(1+O[i]*R*C, 3);


    //x_(k+1) = X_k - alpha*gradient(x_k) => S1
    // double alpha, resistor_old, capacitor_old, S1_C_sum, S1_R_sum;
    // capacitor = capacitor - alpha*S1_C;
    
    // resistor
    // resistor = resistor - alpha*S1_R;
    
    // i = 0;
    // resistor_old, capacitor_old = 0;
    // S1_C_sum, S1_R_sum = 0.0;
    // // alpha = 0.01;
  
    // while (abs(resistor - resistor_old) > 0.001 && abs(capacitor - capacitor_old) > 0.001) {
        
    //     S1_C = 2 * W * (2* pow(omega_array[i], 2) * resistor * capacitor * (resistor - z_real_array[i] * (1 + pow(omega_array[i]*resistor*capacitor, 2))))/pow(1 + pow(omega_array[i]*resistor*capacitor, 2), 3);
    //     // S1_C = 2 * W * ((resistor/1 + pow(resistor*capacitor*omega_array[i], 2))-z_real_array[i])*(2*pow(omega_array[i], 2)*resistor*capacitor)/pow(1+pow(omega_array[i]*resistor*capacitor,2),2);
    //     S1_C_sum = S1_C_sum + S1_C;
        
    //     S1_R = 2 * W * (resistor - z_real_array[i] * (1 + pow(omega_array[i]*resistor*capacitor, 2))) * (-1 + 3 * pow(omega_array[i]*resistor*capacitor, 2))/pow(resistor, 2)*pow(1 + pow(omega_array[i]*resistor*capacitor, 2), 3);
    //     // S1_R = 2*W*((resistor/1+pow(omega_array[i]*resistor*capacitor,2))-z_real_array[i])*(-1+3*pow(omega_array[i]*resistor*capacitor, 2))/pow(resistor,2)*pow(1+pow(omega_array[i]*capacitor*resistor,2),2);
    //     S1_R_sum = S1_R_sum + S1_R;

    //     resistor_old = resistor;
    //     capacitor_old = capacitor;

    //     capacitor = capacitor_old - alpha * S1_C;
    //     resistor = resistor_old - alpha * S1_R;
    //     i++;

    //     cout << "S1_C_sum: " << S1_C_sum << "\n";
    //     cout << "S1_R_sum: " << S1_R_sum << "\n";
    //     cout << "Capacitor: " << capacitor << "\n";
    //     cout << "Resistor: " << resistor << "\n";
        
    // }

    // gradient_descent(1, 0.001, 0.1, 0.1, 10000, omega_array[n_f+1], z_real_array[n_f], z_image_array[n_f])
