import matplotlib.pyplot as plt
import numpy as np




# DADOS
# frequency -> frequência
frequency = np.loadtxt("result.txt", delimiter=" ", usecols=(0))

# data_x -> Z'
z_real = np.loadtxt("result.txt", delimiter=" ", usecols=(1))

# data_y -> Z" 
z_image = np.loadtxt("result.txt", delimiter=" ", usecols=(2))

# TAXA DE APRENDIZADO ( velocidade de descida )
alpha = 0.01


# VALORES INICIAIS PARA c e r
# w0 
c = 0.1

# w1
r = 0.1

def angular_frequency(frequency):
    return 2 * np.pi * frequency

def y_hat(c, r):
    z_real = r * (1+ (angular_frequency(frequency) * c * r)**2)**-1 
    z_image = - (angular_frequency(frequency) * r**2 * c * (1 + (angular_frequency(frequency) * r * c)**2)**-1)
    #plt.plot(z_real,z_image, "bo")
    #plt.savefig("figura.png", dpi=1200)

    return z_real,  z_image


def SC_1(W=1, z_real):
    numerador = (-2*r**3*angular_frequency(frequency)**2*c)
    denominador = (1+ (-angular_frequency(frequency)*r*c)**2)**2 
    sc_1 = 2*W*(y_hat.z_real(c, r) - z_real)*numerador/denominador
    return sc_1


def SC_2(W=1, z_image):
    numerador = -angular_frequency(frequency)*r**2*((-angular_frequency(frequency)*r*c)**2 + 1)
    denominador = (1+ (-angular_frequency(frequency)*r*c)**2)**2
    sc_2 = 2*W*(y_hat.z_image(c, r) - z_image)*numerador/denominador 
    return sc_2


def plot_line(z_real , z_image, c, r):
    x_values = [i for i in range(int(min(z_real)) - 1, int(max(z_real)) + 2)]
    y_values = [y_hat(c, r) for x in x_values]
    plt.plot(x_values,y_values)
    #plt.plot(z_real ,-z_image,'bo')
    #plt.savefig("figura-2.png", dpi=1200)

plot_line(z_real, z_image, c, r) 

def MSE (z_real, z_image, c , r):
    custo = 0
    m = float(len(z_real ))
    for i in range(0, len(z_real )):
        custo += (y_hat(z_real[i], c, r)-z_image[i])**2
    print("custo: ",custo)

        
    return custo/m

#print(MSE(z_real , z_image, c, r))   



def gradient_descent_step(c, r, z_real , z_image, alpha):
    
    erro_c = 0
    erro_r = 0
    W = 1
    m = float(len(z_real ))
    
    for i in range(0,len(z_real )):
        erro_c += 2 * W * y_hat(z_real[i], c, r) - z_image[i]
        erro_r += (y_hat(z_real[i], c, r) - z_image[i]) * z_real [i]
        
    new_c = c - alpha * (1/m) * erro_c
    new_r = r - alpha * (1/m) * erro_r

    return new_c, new_r


epoch = 100

def gradient_descent(c, r, z_real , z_image, alpha, epoch):
    custo = np.zeros(epoch)
    for i in range(epoch):
        c , r = gradient_descent_step(c, r, z_real , z_image, alpha)
        custo[i] = MSE(z_real, z_image , c , r)
        
    return c, r, custo




    
# plt.scatter(data_x, -1*z_image)
# plt.savefig("figura.png", dpi=1200)
