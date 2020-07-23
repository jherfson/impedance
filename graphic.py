import matplotlib.pyplot as plt
import numpy as np

path = "/home/jherfson/Dropbox/Research_Jherfson/Estabilidade termodinâmica e diagramas de fase teóricos/impedancia/python/result.txt"
frequency = np.loadtxt(path, delimiter=" ", usecols=(0))
z_real_txt = np.loadtxt(path, delimiter=" ", usecols=(1))
z_image_txt = np.loadtxt(path, delimiter=" ", usecols=(2))

alpha = 0.01
c = 0.1
r = 0.1

def angular_frequency(frequency):
    return 2 * np.pi * frequency


def y_hat(f, c, r):
    z_real = r * (1+ (angular_frequency(f) * c * r)**2)**-1 
    z_image = - (angular_frequency(f) * r**2 * c * (1 + (angular_frequency(f) * r * c)**2)**-1)
    return z_real, z_image


def plot_line(z_real_txt, z_image_txt, c, r):
    x_values = [i for i in range(int(min(z_real_txt)) - 1, int(max(z_real_txt)) + 2)]
    y_values_z_real = [y_hat(f, c, r) for f in x_values]
    plt.plot(x_values,y_values_z_real, 'r')
    plt.plot(z_real_txt, z_image_txt, 'bo')


def MSE(z_real_txt, z_image_txt):
    custo = 0
    m = float(len(frequency))
    for i in range(0, len(frequency)):
        # custo = custo + (z_real_txt[i] - z_image_txt[i])**2
        custo += (z_real_txt[i] - z_image_txt[i])**2
    print("custo: ",custo)
        
    return custo/m

def gradient_descent_step(c, r, z_real_txt , z_image_txt, alpha):
    
    erro_c = 0
    erro_r = 0
    m = float(len(frequency))
    for i in range(0,len(z_real_txt)):
        W=1
        
        denominador = (1 + (angular_frequency(frequency[i]) * r * c)**2)**2 
        # capacitor
        numerador_c = (-2 * r**3 * angular_frequency(frequency[i])**2 * c)
        sc_1 = 2 * W * (y_hat(frequency[i] ,c, r)[0] - z_real_txt[i]) * numerador_c/denominador
        sc_2 = 2 * W * (y_hat(frequency[i], c, r)[1] - z_image_txt[i]) * (-angular_frequency(frequency[i]) * r**2*((-angular_frequency(frequency[i])* r * c)**2 + 1))/denominador
        
        # Resitor
        numerador_r = (-(angular_frequency(frequency[i]) * c * r)**2 + 1)
        sr_1 = 2 * W * (y_hat(frequency[i], c, r)[0] - z_real_txt[i]) * numerador_r/denominador
        sr_2 = 2 * W * (y_hat(frequency[i], c, r)[1] - z_image_txt[i]) * (-2 * angular_frequency(frequency[i]) * r * c)/denominador 


        erro_c += sc_1 + sc_2
        erro_r += sr_1 + sr_2
        
    new_c = c - alpha * (1/m) * erro_c
    new_r = r - alpha * (1/m) * erro_r

    return new_c, new_r

epoch = 30

def gradient_descent(c, r, z_real_txt, z_image_txt, alpha, epoch):
    custo = np.zeros(epoch)
    for i in range(epoch):
        c , r = gradient_descent_step(c, r, z_real_txt, z_image_txt, alpha)
        custo[i] = MSE(z_real_txt, z_image_txt)
        
    return c, r, custo



print(gradient_descent(c, r, z_real_txt, z_real_txt, alpha, epoch))
