from pylab import *
from numpy import *
import h5py  as h5

def deriv_x(u):
    
    #Get the size of the box
    N=u.shape[0]

    #Define wavenumbers in Fourier space
    k=r_[0:N//2,-N//2:0]
    kx=k.reshape((-1, 1, 1))
    ky=k.reshape(( 1,-1, 1))
    kz=k[0:(N//2+1)].reshape(( 1, 1,-1))
    
    du=irfftn(1j*kx*rfftn(u))
    
    return du

def deriv_y(u):
    
    N=u.shape[0]

    #k=fftshift(arange(0,N)-N//2)
    k=r_[0:N//2,-N//2:0]
    kx=k.reshape((-1, 1, 1))
    ky=k.reshape(( 1,-1, 1))
    kz=k[0:(N//2+1)].reshape(( 1, 1,-1))
    
    #Compute the derivative using FFT 
    du=irfftn(1j*ky*rfftn(u))
    
    return du
  
def deriv_z(u):
    
    N=u.shape[0]

    #k=fftshift(arange(0,N)-N//2)
    k=r_[0:N//2,-N//2:0]
    kx=k.reshape((-1, 1, 1))
    ky=k.reshape(( 1,-1, 1))
    kz=k[0:(N//2+1)].reshape(( 1, 1,-1))
    
    #Compute the derivative using FFT 
    du=irfftn(1j*kz*rfftn(u))
    
    #Return
    return du

def deriv(u,flag):

    #Get the size of the box
    N=u.shape[0]

    #Define wavenumbers in Fourier space
    k=r_[0:N//2,-N//2:0]
    kx=k.reshape((-1, 1, 1))
    ky=k.reshape(( 1,-1, 1))
    kz=k[0:(N//2+1)].reshape(( 1, 1,-1))
    
    #Choose direction
    if(flag==0):
        kk=kx
    if(flag==1):
        kk=ky
    if(flag==2):
        kk=kz
        
    #Compute the derivative 
    dux=irfftn(1j*kk*rfftn(u))
    
    #Return
    return dux
  
def laplacian(u):

    #Get the size of the box
    N=u.shape[0]

    #Define wavenumbers in Fourier space
    k=r_[0:N//2,-N//2:0]
    kx=k.reshape((-1, 1, 1))
    ky=k.reshape(( 1,-1, 1))
    kz=k[0:(N//2+1)].reshape(( 1, 1,-1))
    
    #Laplace operator
    kk=kx**2 + ky**2 + kz**2
    kk[0,0,0]=1
    
    #Invert the Laplace operator 
    lap=irfftn(-1.0/kk*rfftn(u))
    
    #Return
    return lap

def randomise(u):
    
    #Get size of the box 
    N=u.shape[0]
    
    #Generate random numbers in 0-2*pi
    phi=2*pi*rand(N,N,N//2+1)
    rand_unit=exp(1j*phi)
    
    #Shift Fourier modes randomly
    u_rand=irfftn(rand_unit*rfftn(u))
    
    #Return
    return u_rand
    
