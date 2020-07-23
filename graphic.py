import matplotlib.pyplot as plt
import numpy as np

path = "/home/jherfson/Dropbox/Research_Jherfson/Estabilidade termodinâmica e diagramas de fase teóricos/impedancia/python/result.txt"
frequency = np.loadtxt(path, delimiter=" ", usecols=(0))
z_real_txt = np.loadtxt(path, delimiter=" ", usecols=(1))
z_image_txt = np.loadtxt(path, delimiter=" ", usecols=(2))

alpha = 0.01
c = 0.1
r = 0.1
max_itera = 100
step_tolerance = 0.001

def angular_frequency(frequency):
    return 2 * np.pi * frequency


def y_hat(frequency, c, r):
    z_real = r * (1+ (angular_frequency(frequency) * c * r)**2)**-1 
    z_image = - (angular_frequency(frequency) * r**2 * c * (1 + (angular_frequency(frequency) * r * c)**2)**-1)
    return z_real, z_image


def gradient_descent(c, r, z_real_txt ,z_image_txt, alpha):
    
    step_c = 0
    step_r = 0
    for i in range(max_itera):
        W=1
        print(f'Frequência: {frequency[i]}')

        denominador = (1 + (angular_frequency(frequency[i]) * r * c)**2)**2 
        # capacitor
        numerador_c = (-2 * r**3 * angular_frequency(frequency[i])**2 * c)
        sc_1 = 2 * W * (y_hat(frequency[i] ,c, r)[0] - z_real_txt[i]) * numerador_c/denominador
        sc_2 = 2 * W * (y_hat(frequency[i], c, r)[1] - z_image_txt[i]) * (-angular_frequency(frequency[i]) * r**2*((-angular_frequency(frequency[i])* r * c)**2 + 1))/denominador
        sc = sc_1 + sc_2
        # Resitor
        numerador_r = (-(angular_frequency(frequency[i]) * c * r)**2 + 1)
        sr_1 = 2 * W * (y_hat(frequency[i], c, r)[0] - z_real_txt[i]) * numerador_r/denominador
        sr_2 = 2 * W * (y_hat(frequency[i], c, r)[1] - z_image_txt[i]) * (-2 * angular_frequency(frequency[i]) * r * c)/denominador 
        sr = sr_1 + sr_2

        step_c = alpha * sc        
        c -= step_c
        print(f'Novo C:{c}')

        step_r = alpha * sr
        r -= step_r
        print(f'Novo R:{r}')

        if abs((step_c <= step_tolerance and step_r <= step_tolerance)):
            break

    return c, r

# epoch = 30

# def gradient_descent(c, r, z_real_txt, z_image_txt, alpha, epoch):
#     custo = np.zeros(epoch)
#     for i in range(epoch):
#         c , r = gradient_descent_step(c, r, z_real_txt, z_image_txt, alpha)
#         custo[i] = MSE(z_real_txt, z_image_txt)
        
#     return c, r, custo



print(gradient_descent(c, r, z_real_txt, z_image_txt, alpha))
