# import matplotlib.pyplot as plt
import numpy as np
import os
path = os.path.join(os.path.abspath('.'), 'result.txt')
impedancia = np.loadtxt(path, delimiter=" ")


frequency = [impedancia[i][0] for i in range(len(impedancia))]
z_real_txt = [impedancia[i][1] for i in range(len(impedancia))]
z_imag_txt = [impedancia[i][2] for i in range(len(impedancia))]

f = np.array(frequency)
#z_exp = 1/(1/5e4+1/(1/(complex(0, 1)*2*np.pi*f*5e-8)))
#z_re_exp = z_exp.real
#z_im_exp = z_exp.imag
z_re_exp = np.array(z_real_txt)
z_im_exp = np.array(z_imag_txt)


#alpha = 0.001
C = 4e-7  # capacitor -> fara
R = 100000  # resitor -> Ohm

def gradient_descent(C, R, f, z_re_exp, z_im_exp):
    i = complex(0, 1)
    w = 2*np.pi*f
    Ro = R
    Co = C
    r = R/Ro
    c = C/Co


    W = 1
    C_old = 2 * C
    R_old = 2 * R
    delta_r = 1
    delta_c = 1
    controller = ''
    count = 0
    outputfile = open("test.txt", "w")
    while (abs(R_old - R)/R_old > 0.0001 or abs(C_old - C)/C_old > 0.0001) and count < 1000:
        #print(f'cycle begining; R={R},C={C}')
        z_c = 1/(i * w * C)
        
        z_r = R
        z = 1/(1/z_r + 1/z_c)
        z_re = z.real
        z_im = z.imag

        tau = R * C
        factor1 = 1 - (w * tau)**2
        factor2 = -2 * w * tau
        factor3 = w * R**2
        denominador = (1 + (w * tau)**2)**2
        s1 = np.sum(W*(z_re-z_re_exp)**2)
        s2 = np.sum(W*(z_im-z_im_exp)**2)
        s = s1+s2
        #print(f'sum of squares: {s}')

        # capacitor
        # print(f'numerador_c:{numerador_c}')
        sc_1 = 2 * W * (z_re - z_re_exp) * Co*factor2*factor3/denominador
        sc_2 = 2 * W * (z_im - z_im_exp) * Co*(-factor1 * factor3)/denominador
        sc = sc_1 + sc_2
        #print(f'somas (sc):{np.sum(sc_1)},{np.sum(sc_2)}')
        sc_sum = np.sum(sc)

        if sc_sum == 0:
            alpha_c = 1
        else:
            ratio_c = min(0.90*delta_c*c,0.5*c) if controller == 'c' else 0.5*c 
            alpha_c = min(1, ratio_c/abs(sc_sum))

        # Resitor

        sr_1 = 2 * W * (z_re - z_re_exp) * Ro*factor1/denominador
        sr_2 = 2 * W * (z_im - z_im_exp) * Ro*factor2/denominador
        sr = sr_1 + sr_2
        #print(f'somas (sr):{np.sum(sr_1)},{np.sum(sr_2)}')
        sr_sum = np.sum(sr)

        if sr_sum == 0:
            alpha_r = 1
        else:
            ratio_r = min(0.90*delta_r*r,0.5*r) if controller == 'r' else 0.5*r 
            alpha_r = min(1, ratio_r/abs(sr_sum))

        alpha = min(alpha_c, alpha_r)
        if alpha_c < alpha_r:
            controller = 'c'
        else:
            controller = 'r'


        
        print(f'Alpha: {alpha},alpha_r:{alpha_r},alpha_c:{alpha_c}')
        
        c_old = c
        C_old = c_old*Co
        print(f'c: {c}; sc_sum: {sc_sum}')
        c -= alpha * sc_sum
        delta_c = abs(c-c_old)
        C = c*Co
        print(f'Novo c:{c}, Novo C:{C}')
        
        r_old = r
        R_old = r_old*Ro
        print(f'r: {r}; sr_sum: {sr_sum}')
        r -= alpha * sr_sum
        delta_r = abs(r-r_old)
        R = r*Ro 
        print(f'Novo r:{r}, Novo R:{R}\n')
        
        Co = C 
        Ro = R
        
        c = 1
        r = 1
        
        # outputfile.write(" \n%s" % str(alpha)+';   '+str(C)+';   '+str(R)+';   '+str(sc_sum)+';   '+str(sr_sum))

        count += 1

    print(f'Contador:{count}, C_old:{C_old}, R_old:{R_old}\n')
    print(r'δC: ', C-C_old)
    print(r'δR: ', R-R_old)
    return C, R, sc_sum, sr_sum


print(gradient_descent(C, R, f, z_re_exp, z_im_exp))
