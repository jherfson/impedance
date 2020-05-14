// Implementação do sistema pra impedância com o método de gradiente
// versão: 0.1

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include <fstream>


using std::complex;
// using std::cout;
// using std::endl;
// using std::cin;

using namespace std; 



//using COMPLX = complex<double>;

// protótipo da função da frequência angular
double angular_frequency(double frequency); // protótipo da função

// impedância resistor
double z_resistor(double R);

// impedância capacitor
double z_capacitor(double capacitor);

// impedância CPE
double z_CPE(double frequency, double T);

void output(string filename);

// função principal
int main() 
{
    double f, T, r;
    string filename;
    // cout << "Digite o valor frequência:\n";
    // cin >> f;
    // cout << "Digite o Periodo:\n";
    // cin >> T; 
    // cout << "Digite o valor do resistor:\n";
    // cin >> r; 
    // cout << "Valor da frequência angular " << angular_frequency(f) << "\n";
    // cout << "Valor da frequência: "<< f << "\n";
    
    // cout << "-----------------------------------\n";

    // cout << "Valor z_CPE " << z_CPE(f, T) << "\n";
    // cout << "Valor do Período: " << T << endl;

    cout << "Digite o nome do Arquivo: ";
    cin >> filename;
    output(filename);



    // cout << "-----------------------------------\n";

    // cout << "Valor Resistor " << z_resistor(r) << "\n";
    // cout << "Valor do Período: " << r << endl;


    return 0;

}

// função da frequência angular
double angular_frequency(double frequency) 
{
    double w;
    w = 2*M_PI*frequency;
    return w;

}

// funçõ da impedância resistor
// double z_resistor(double resistor)
// {
//     //double z;
//     complex<double> z(resistor, 0);
//     return 0;
     
// }

// função da impedância capacitor
double z_capacitor(double capacitor, double w)
{
    double z;
    complex<double> i(0, 1);
    //z = 1/abs(i)*angular_frequency(w)*capacitor;
    z = pow(abs(i)*angular_frequency(w)*capacitor, -1);


    return z;
}

// função da impedância CPE
double z_CPE(double frequency, double T)
{
    double w, z;
    int p = 1;
    w = angular_frequency(frequency);
    complex<double> i(0, 1);
    // z = 1/T*pow((abs(i)*w), p);
    z = pow(T*pow((abs(i)*w), p), -1); 

    return z;

}

// Leitura do arquivo com os dados da impedância
void output(string filename)
{
    ifstream txtFile;
    string STRING;


    int tam = 0;
    string f[tam], z_real[tam], z_im[tam];

    txtFile.open(filename);

    while(!txtFile.eof())
    {
        tam++;
        //getline(txtFile, STRING);
        getline(txtFile, STRING,' ');

        f[tam] = STRING;
        z_real[tam] = STRING;
        z_real[tam] = STRING;

        




        cout << f[tam] ;
        cout << z_real[tam] ;
        cout << z_real[tam] << endl;

    }

    txtFile.close();

}

double sum_erro_impedance(double f, double z_real, double z_img)
{

}




// void output(string filename)
// {
//     FILE *file;
//     //Quantidade de linhas (pode servir pra algum calculo estatistico)
//     int quant=0;
//     float x[quant], y[quant], z[quant];
//     //Tamanho de caracteres encontrados em uma linha (chute)
//     char linha[50];

//     //Não foi feito aqui aquele teste pra saber se o arquivo foi lido corretamente (não precisei)
//     file = fopen("dados.txt", "r");
    
//     //Enquanto não alcançar o fim do arquivo, faça o seguinte linha por linha:
//     while(fgets(linha, sizeof(linha), file)) 
//     {
//         x[quant] = atof(strtok())
//         x[quant] = atof(strtok(linha, " "));
//         y[quant] = atof(strtok(NULL, " "));
//         quant++;    
//     }
//     fclose(file);
    
//     system("pause");
    
//     return 0;
// }





