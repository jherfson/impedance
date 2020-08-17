# import matplotlib.pyplot as plt
import os
import numpy as np
from scipy.optimize import line_search
import matplotlib.pyplot as plt

path = os.path.join(os.path.abspath('.'), 'result.txt')
impedancia = np.loadtxt(path, delimiter=" ")


frequency = [impedancia[i][0] for i in range(len(impedancia))]
z_real_txt = [impedancia[i][1] for i in range(len(impedancia))]
z_imag_txt = [impedancia[i][2] for i in range(len(impedancia))]

f = np.array(frequency)
z_exp = 1/(1/5e4+1/(1/(complex(0, 1)*2*np.pi*f*5e-8)))
z_re_exp = z_exp.real
z_im_exp = z_exp.imag
# z_re_exp = np.array(z_real_txt)
# z_im_exp = np.array(z_imag_txt)


#alpha = 0.001
C = 4e-7  # capacitor -> fara
R = 100000  # resitor -> Ohm

# def angular_frequency(f):
#     return 2 * np.pi * f


def sum_errors(z_re_exp, z_im_exp, f):
    i = complex(0, 1)
    w = 2 * np.pi * f
    z_c = 1/(i * w * C)
    z_r = R
    z = 1/(1/z_r + 1/z_c)

    return z.real, z.imag

def gradient():
    

print(sum_errors(z_re_exp, z_im_exp, f))

