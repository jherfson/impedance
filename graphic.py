import matplotlib.pyplot as plt
import numpy as np


# DADOS
# data_x -> Z'
data_x= np.loadtxt("result.txt", delimiter=" ", usecols=(1))

# data_y -> Z" 
data_y = np.loadtxt("result.txt", delimiter=" ", usecols=(2))

# TAXA DE APRENDIZADO ( velocidade de descida )
alpha = 0.01


# VALORES INICIAIS PARA c e r
# w0 
c = 0.1

# w1
r = 0.1

def y_hat(x, c, r):
    return c + r * data_x

print(y_hat(1.5, c, r))


def plot_line(data_x, data_y, c, r):
    x_values = [i for i in range(int(min(data_x)) - 1, int(max(data_x)) + 2)]
    y_values = [y_hat(x, c, r) for x in x_values]
    plt.plot(x_values,y_values,'r')
    plt.plot(data_x,data_y,'bo')
    #plt.savefig("figura-2.png", dpi=1200)

plot_line(data_x, data_y, c, r) 

def MSE (data_x, data_y , c , r):
    custo = 0
    m = float(len(data_x))
    for i in range(0, len(data_x)):
        custo += (y_hat(data_x[i], c, r)-data_y[i])**2
        
    return custo/m

#print(MSE(data_x, data_y, c, r))   


def gradient_descent_step(c, r, data_x, data_y, alpha):
    
    erro_c = 0
    erro_r = 0
    m = float(len(data_x))
    
    for i in range(0,len(data_x)):
        erro_c += y_hat(data_x[i], c, r) - data_y[i]
        erro_r += (y_hat(data_x[i], c, r) - data_y[i]) * data_x[i]
        
    new_c = c - alpha * (1/m) * erro_c
    new_r = r - alpha * (1/m) * erro_r

    return new_c, new_r


epoch = 100

def gradient_descent(c, r, data_x, data_y, alpha, epoch):
    custo = np.zeros(epoch)
    for i in range(epoch):
        c , r = gradient_descent_step(c, r, data_x, data_y, alpha)
        custo[i] = MSE(data_x, data_y , c , r)
        
    return c, r, custo


print(gradient_descent(c, r, data_x, data_y, alpha, epoch))

    
# plt.scatter(data_x, -1*data_y)
# plt.savefig("figura.png", dpi=1200)