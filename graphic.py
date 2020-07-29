# import matplotlib.pyplot as plt
import numpy as np
import os

path = os.path.join(os.path.abspath('.'), 'result.txt')
impedancia = np.loadtxt(path, delimiter=" ")


frequency = [impedancia[i][0] for i in range(len(impedancia))]
z_real_txt = [impedancia[i][1] for i in range(len(impedancia))]
z_imag_txt = [impedancia[i][2] for i in range(len(impedancia))]

f = np.array(frequency)
z_exp = 1/(1/1e4+1/(1/(complex(0, 1)*2*np.pi*f*1e-8)))
z_re_exp = z_exp.real
z_im_exp = z_exp.imag
# z_re_exp = np.array(z_real_txt)
# z_im_exp = np.array(z_imag_txt)


#alpha = 0.001
Co = 2e-8  # capacitor -> fara
Ro = 20000  # resitor -> Ohm

# c = 0.001
# r = 1
# z_c = 1/(i*w*c)
# z_r = r
# z = 1/( 1/z_r + 1/z_c )
# z_re = z.real
# z_im = z.imag

# tau = r*c


# def angular_frequency(frequency):
#     return 2 * np.pi * frequency


# def impedance(frequency, c, r):
#     z_real = r * (1 + (angular_frequency(frequency) * c * r)**2)**-1
#     z_imag = - (angular_frequency(frequency) * r**2 * c *
#                 (1 + (angular_frequency(frequency) * r * c)**2)**-1)
#     return complex(z_real, z_imag)
# Ro -> chute inicial R => r = R/Ro
# Co -> chute incial C  => c = C/Co
R = Ro
C = Co

r = R/Ro
c = C/Co


def gradient_descent(c, r, f, z_re_exp, z_im_exp):
    i = complex(0, 1)
    w = 2*np.pi*f
    C = c
    R = r


    W = 1
    c_old = 2 * c
    r_old = 2 * r
    count = 0

    while (abs(r_old - r)/r_old > 0.001 or abs(c_old - c)/c_old > 0.001) and count < 10:
        # z_c = 1/(i * w * c)
        z_c = 1/(i * w * C)
        # z_r = r
        z_r = R
        z = 1/(1/z_r + 1/z_c)
        z_re = z.real
        z_im = z.imag

        tau = r * c
        factor1 = 1 - (w * tau)**2
        factor2 = -2 * w * tau
        factor3 = w * r**2
        denominador = (1 + (w * tau)**2)**2
        s1 = np.sum(W*(z_re-z_re_exp)**2)
        s2 = np.sum(W*(z_im-z_im_exp)**2)
        s = s1+s2
        #print(f'sum of squares: {s}')

        # capacitor
        # print(f'numerador_c:{numerador_c}')
        sc_1 = 2 * W * (z_re - z_re_exp) * factor2*factor3/denominador
        sc_2 = 2 * W * (z_im - z_im_exp) * (-factor1 * factor3)/denominador
        sc = sc_1 + sc_2
        # print(f'somas:{np.sum(sc_1)},{np.sum(sc_2)}')
        sc_sum = np.sum(sc)

        if sc_sum == 0:
            alpha_c = 1
        else:
            alpha_c = min(1, 0.5*c/abs(sc_sum))

        # Resitor

        sr_1 = 2 * W * (z_re - z_re_exp) * factor1/denominador
        sr_2 = 2 * W * (z_im - z_im_exp) * factor2/denominador
        sr = sr_1 + sr_2
        sr_sum = np.sum(sr)

        if sr_sum == 0:
            alpha_r = 1
        else:
            alpha_r = min(1, 0.5*r/abs(sr_sum))

        alpha = min(alpha_c, alpha_r)
        # alpha = 0.01
        print(f'Alpha: {alpha}')
        # step_c = alpha * sc
        # c -= step_c
        c_old = c
        print(f'c: {c}; sc_sum: {sc_sum}')
        c -= alpha * sc_sum
        # print(f'Step C:{step_c}')
        print(f'Novo c:{c}')
        # step_r = alpha * sr
        # r -= step_r
        r_old = r
        print(f'r: {r}; sr_sum: {sr_sum}')
        r -= alpha * sr_sum
        C = c*Co
        R = r*Ro
        print(f'Novo r:{r}')
        print(f'C = rCo: {c*Co}')
        print(f'R = rRo: {r*Ro}\n')
        C = c*Co
        R = c*Co
        count += 1

    print(f'Contador:{count}, C_old:{c_old}, R_old:{r_old}\n')
    print(r'$\Delta$c: ', (c_old-c))
    print(r'$\Delta$r: ', (r_old-r))
    return c, r, sc_sum, sr_sum


print(gradient_descent(c, r, f, z_re_exp, z_im_exp))
