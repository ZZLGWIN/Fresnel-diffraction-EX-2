# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 15:05:45 2018

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

#function numerically intergrates a function using the simpsons rule 
def intergrate(func, N , upperbound, lowerbound):
    
    Sum=func(upperbound)+func(lowerbound)
    delta_x=(upperbound-lowerbound)/N
   
    #N has to be even 
    if N%2!=0:
        return print('N must be even')
   #there are 2 series an even series and an odd series
   #the diffrent series have to be multiplyed by diffrent coefficants
    else:    
        n=1
        #odd series
        while n<N:
            n+=2
            Sum+=4*func(lowerbound+n*delta_x)
        n=2
        #even series
        while n<N-1:
            n+=2
            Sum+=2*func(lowerbound+n*delta_x)
        return Sum*delta_x/3
    
def Diffraction(wavelength,ScreenSeperation,Numpoints,xPrimeMin,xPrimeMax,yPrimeMin,yPrimeMax):
    k=2*np.pi/(wavelength*1e-9)
    z=ScreenSeperation*(1e-3)
    #Screen size
    xymin=-1e-4
    xymax=1e-4
    step=(xymax-xymin)/(Numpoints-1)
    x_val=np.zeros(Numpoints)
    y_val=np.zeros(Numpoints)
    X=np.zeros(Numpoints)
    intensity=np.zeros((Numpoints,Numpoints))
    E0=8.854187817e-12


    for i in range(Numpoints):
        x_val[i]=xymin+(i*step)
        X=intergrate(lambda x:np.exp((1j*k)/(2*z)*(x_val[i]-x)**2),1000,xPrimeMax,xPrimeMin)
        for j in range(Numpoints):
            y_val[j]=xymax-(j*step)
            temp = (k/(2*np.pi*z))*X*intergrate(lambda y1:np.exp(1j*k/(2*z)*(y_val[j]-y1)**2),50,yPrimeMax,yPrimeMin)
            intensity[i,j]=(3e8)*(E0)*(temp*np.conjugate(temp))
    
    return intensity

#The following creats a loop that allows the user to choose which part of the exercise they want to look at.
MyInput = '0'
while MyInput != 'q':
    MyInput=input('Enter a choice "a","b","c" or "q" to quit: ')
    if MyInput=='a':
        print('You have choosen part a.')
        print('As a test of my simpsons function the simpsons intergral will be performed on sin(x) between 0 and pi.')
        print(intergrate(lambda x:np.sin(x),1000000,np.pi,0))
  
    elif MyInput=='b':
        def slit_diffraction(Wavelength=500,ScreenSeperation=1,Length=2,xmin=2,xmax=2,plot=False):
            
            Length=Length*(1e-3)
            k=2*np.pi/(Wavelength*(1e-9))
            z=ScreenSeperation*(1e-3)
            x1max=Length/2
            x1min=-Length/2
            
            x_domain=[]
            X=[]
            
            for i in range(int(xmin*100),int(xmax*100),1):
                  x_domain+=[i/100]
                  temp=intergrate(lambda x:np.exp((1j*(k/2*z))*((i/100)-x)**2),10000,x1max,x1min)
                  X+=[temp*np.conjugate(temp)]
                
            if plot==False:
                return X
            if plot==True:
                plt.title('single slit diffraction')
                plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
                plt.legend('screen seperation=' + str(z))
                plt.ylabel('Intensity (W/m^2)')
                plt.xlabel('Position (m)')
                plt.plot(x_domain,X)
                plt.show()
                print('Diffraction pattern for wavelength=' + str(Wavelength) + 'nm, screen seperation=' +str(ScreenSeperation) + 'mm and Slit length=' +str(Length*(1e3)) +'(mm).')
                
                
            
            
        print('You have choosen part b. \nIn this section the diffraction pattern is generated for single slit diffraction.')
        W=float(input('Choose a wavelength (nm): ')) 
        #500nm is a good choice
        Z=float(input('select a screen seperation (mm): '))
        #5mm
        L=float(input('select a appeture size (mm): '))
        #1mm
        Min=float(input('select the lowerbound of the x domain: '))
        #-1
        Max=float(input('select the upperbound of the x domain: '))
        #1
        print('Please wait a moment while the intensity plots are generated. \n')
        slit_diffraction(Wavelength=W,ScreenSeperation=Z,Length=L,xmin=Min,xmax=Max,plot=True)
    
    elif MyInput=='c':
        print('You have choosen part c.')
        wavelength=float(input('Choose a wavelength (nm): '))
        seperation=float(input('select a screen seperation (mm): '))
        print('Please wait a moment while the intensity plots are generated. \n')
        
        Intensity=Diffraction(wavelength,seperation,500,-1e-4,1e-4,-1e-4,1e-4)
        
        plt.imshow(Intensity,cmap='inferno')
        plt.colorbar()
        plt.show()
    
    elif MyInput != 'q':
        print('Input not reconised please try again')
        
        
print('Program ended \nGoodbye')
          





