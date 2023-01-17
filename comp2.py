import numpy as np
import math
from math import pi
import pandas as pd



tv_f = -0.098363164308346596735
tv_g = 0.11033333333333333333
tv_h = 1.7688741449316336902  # True value of the integral by Wolfram Alpha

f = lambda x: np.cos(pi*x)

def g(x):
    if x>=0:
        y=x**2
    else:
        y=-1*x**2
        
    return y

h = lambda x: np.exp(-x^2/2)

def methodf(n):
    a = -1.0
    b = 1.1
    h = float(b - a) / n
    f = lambda x: np.cos(pi*x)
    k = 0.0
    
    for i in range (1,n):
        k= k+ h*f(a+i*h)
    
    v_CTR = k+(h/2)*(f(a)+f(b))
    
    
    k = 0.0
    x=a + h
    if (n % 2) == 0:
        for i in range(1,(n//2)+1):
            k += 4*f(x)
            x += 2*h

        x = a + 2*h
        for i in range(1,n//2):
            k += 2*f(x)
            x += 2*h
        
        
    v_CSR = (h/3)*(f(a)+f(b)+k)
        
    
    
    
    # By applying the change of Interval
    # Gaussian
    v_Gau = 0.0
    
    x0=-np.sqrt(3/5)
    x1=0.0
    x2=np.sqrt(3/5)

    for i in range(n):
        upper = a+ (i+1)*h
        lower = a+ i*h
        
        sub = (upper-lower)/2
        add = (upper+lower)/2
        
        v_Gau=v_Gau+ (5/9*f(x0*sub + add)+8/9* f(x1*sub + add) + 5/9 * f(x2 * sub + add))
        
        
    v_Gau = h/2 * v_Gau
                
    #v_Gau = 0.05*(5/9*f(x0*0.05+1.05)+8/9*f(x1*0.05+1.05)+5/9*(x2*0.05+1.05))
    
    
    
    
    e_CTR = abs(tv_f-v_CTR)
    e_CSR = abs(tv_f-v_CSR)
    e_Gau = abs(tv_f-v_Gau)
    
    return e_CTR, e_CSR, e_Gau

def methodg(n):
    a = -1.0
    b = 1.1
    h = float(b - a) / n 
    k = 0.0
    
    for i in range (1,n):
        k= k+ h*g(a+i*h)
    
    v_CTR = k+(h/2)*(g(a)+g(b))
    
    
    k = 0.0
    x= a + h
    if (n % 2) == 0:
        for i in range(1,(n//2)+1):
            k += 4*g(x)
            x += 2*h

        x = a + 2*h
        for i in range(1,n//2):
            k += 2*g(x)
            x += 2*h
        
        
    v_CSR = (h/3)*(g(a)+g(b)+k)
        
        
    
    # By applying the change of Interval
    # Gaussian
    v_Gau = 0.0
    
    x0=-np.sqrt(3/5)
    x1=0.0
    x2=np.sqrt(3/5)

    for i in range(n):
        upper = a+ (i+1)*h
        lower = a+ i*h
        
        sub = (upper-lower)/2
        add = (upper+lower)/2
        
        v_Gau=v_Gau+ (5/9*g(x0*sub + add)+8/9* g(x1*sub + add) + 5/9 * g(x2 * sub + add))
        
        
    v_Gau = h/2 * v_Gau
        
    e_CTR = abs(tv_g-v_CTR)
    e_CSR = abs(tv_g-v_CSR)
    e_Gau = abs(tv_g-v_Gau)
    
    return e_CTR, e_CSR, e_Gau


def methodh(n):
    a = -1.0
    b = 1.1
    h = float(b - a) / n
    hf = lambda x: np.exp((-x**2)/2)
    k = 0.0
    
    for i in range (1,n):
        k= k+ h*hf(a+i*h)
    
    v_CTR = k+(h/2)*(hf(a)+hf(b))
    
    
    k = 0.0
    x=a + h
    if (n % 2 == 0):
        for i in range(1,n//2 + 1):
            k += 4*hf(x)
            x += 2*h

        x = a + 2*h
        for i in range(1,n//2):
            k += 2*hf(x)
            x += 2*h
        
        
        
    v_CSR = (h/3)*(hf(a)+hf(b)+k)
    
    # By applying the change of Interval
    # Gaussian
    v_Gau = 0.0
    
    x0=-np.sqrt(3/5)
    x1=0.0
    x2=np.sqrt(3/5)

    for i in range(n):
        upper = a+ (i+1)*h
        lower = a+ i*h
        
        sub = (upper-lower)/2
        add = (upper+lower)/2
        
        v_Gau=v_Gau+ (5/9*hf(x0*sub + add)+8/9* hf(x1*sub + add) + 5/9 * hf(x2 * sub + add))
        
        
    v_Gau = h/2 * v_Gau
    
    e_CTR = abs(tv_h-v_CTR)
    e_CSR = abs(tv_h-v_CSR)
    e_Gau = abs(tv_h-v_Gau)
    
    return e_CTR, e_CSR, e_Gau




def run():
    df = pd.DataFrame (0, columns = [["Trapezoidal","Trapezoidal","Trapezoidal","Simpson","Simpson","Simpson","Gauss3","Gauss3","Gauss3"],
                                    ["cosine","C1 function","normal dist","cosine","C1 function","normal dist","cosine","C1 function","normal dist"]],
                          index = [2**0,2**1,2**2,2**3,2**4,2**5,2**6,2**7,2**8,2**9]) 
    ecos = np.zeros([10,3])
    ec1 = np.zeros([10,3])
    endist = np.zeros([10,3])
    
    for i in range(0,10):
        
        n = 2**i
        
        x,y,z = methodf(n)
        ecos[i,0] = x
        ecos[i,1] = y
        ecos[i,2] = z
        
        x,y,z = methodg(n)
        ec1[i,0] = x
        ec1[i,1] = y
        ec1[i,2] = z
        
        x,y,z = methodh(n)
        endist[i,0] = x
        endist[i,1] = y
        endist[i,2] = z
        
    for  i in range (1,10):
        df.iloc[i,0] = (np.log(ecos[i-1,0]/ecos[i,0]))/(np.log(2))
        df.iloc[i,1] = (np.log(ec1[i-1,0]/ec1[i,0]))/(np.log(2))
        df.iloc[i,2] = (np.log(endist[i-1,0]/endist[i,0]))/(np.log(2))
        df.iloc[i,3] = (np.log(ecos[i-1,1]/ecos[i,1]))/(np.log(2))
        df.iloc[i,4] = (np.log(ec1[i-1,1]/ec1[i,1]))/(np.log(2))
        df.iloc[i,5] = (np.log(endist[i-1,1]/endist[i,1]))/(np.log(2))
        df.iloc[i,6] = (np.log(ecos[i-1,2]/ecos[i,2]))/(np.log(2))
        df.iloc[i,7] = (np.log(ec1[i-1,2]/ec1[i,2]))/(np.log(2))
        df.iloc[i,8] = (np.log(endist[i-1,2]/endist[i,2]))/(np.log(2))
        
        
    df_2 = pd.DataFrame (0, columns = [["cosine","cosine","cosine","C1 function","C1 function","C1 function","normal dist","normal dist","normal dist"],
                                    ["Trapezoidal","Simpson","Gauss3","Trapezoidal","Simpson","Gauss3","Trapezoidal","Simpson","Gauss3"]],
                          index = [2**0,2**1,2**2,2**3,2**4,2**5,2**6,2**7,2**8,2**9]) 
    for  i in range (1,10):
        df_2.iloc[i,0] = df.iloc[i,0]
        df_2.iloc[i,1] = df.iloc[i,3]
        df_2.iloc[i,2] = df.iloc[i,6]
        df_2.iloc[i,3] = df.iloc[i,1]
        df_2.iloc[i,4] = df.iloc[i,4]
        df_2.iloc[i,5] = df.iloc[i,7]
        df_2.iloc[i,6] = df.iloc[i,2]
        df_2.iloc[i,7] = df.iloc[i,5]
        df_2.iloc[i,8] = df.iloc[i,8]
    
    
    return df,df_2