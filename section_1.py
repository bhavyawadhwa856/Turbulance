from pylab import *
import h5py  as h5
from utils import *

#Use matplotlib interactively
ion()

#Read directory
path_r="./"

#Read file
filename="field_1.h5"

#Open the file
f = h5.File(path_r+filename, 'r')

#Read the three velocity components
u = f['u'][:]
v = f['v'][:]
w = f['w'][:]

#Read the kinematic viscosity
l=2.0/(sqrt(2.0)/3.0*128)
nu=(l**2/l**(2.0/3.0))

#Plot a pencil
figure(1)
#Get the size of the box
N=u.shape[0]
#Get the mesh: the size of the box is 2*pi
x=2*pi*r_[0:N]/N
plot(x,u[1,1,:])
xlabel('x')
ylabel('u')

#Plot a plane
figure(2)
#Get the mesh: the size of the box is 2*pi
xx,yy=meshgrid(x,x)
#Plot a cut of the flow
pcolor(xx,yy,u[1,:,:])
xlabel('x')
ylabel('y')
axis("equal")

###################################################
##ASSIGNMENT 1: What about the small scales? #######
###################################################
##Calc x derivative using deriv_x() and plot them

dux=deriv_x(u)

figure(3)
plot(x,dux[1,1,:])

#Plot a plane
figure(4)
pcolor(xx,yy,dux[1,:,:])
xlabel('x')
ylabel('y')
axis("equal")


###################################################
##ASSIGNMENT 2: Is there turbulence structure?######
###################################################
#Randomise the flow using 
#u_rand=randomise(u), plot and compare

u_rand=randomise(u)
dux_rand=deriv_x(u_rand)

figure(5)
pcolor(xx,yy,dux_rand[1,:,:])
xlabel('x')
ylabel('y')
axis("equal")

###################################
##ASSIGNMENT 3: Mass conservation ##
###################################
#Calculate the divergence of the flow: div=d_x u + d_y v + d_z w

U =zeros((3,N,N,N))

U[0,:,:,:]=u
U[1,:,:,:]=v
U[2,:,:,:]=w

DU=zeros((N,N,N))

for i in range(0,3):
    DU[:,:,:]=DU[:,:,:] + deriv(U[i,:,:,:],i)

########################################################
##ASSIGNMENT 4: Calculate the velocity gradient tensor ##
########################################################
#Calculate \partial_i u_j using deriv_iter()

GU=zeros((3,3,N,N,N))

for i in range(0,3):
    for j in range(0,3):
        GU[i,j,:,:,:]=deriv(U[i,:,:,:],j)
        


########################################################
##ASSIGNMENT 5: Calculate the non-linear terms         ##
########################################################

NL=zeros((3,N,N,N))

for i in range(0,3):
    for j in range(0,3):
        NL[i,:,:,:]=NL[i,:,:,:] + U[j,:,:,:]*GU[i,j,:,:,:]


########################################################
##ASSIGNMENT 6: Calculate the pressure and its grad    ##
########################################################
#Use the function laplacian(u)! Plot pressure

DNL=zeros((N,N,N))
P  =zeros((N,N,N))

for i in range(0,3):
    DNL[:,:,:]=DNL[:,:,:] - deriv(NL[i,:,:,:],i)
        
P=laplacian(DNL) 

GP=zeros((3,N,N,N))

for i in range(0,3):
    GP[i,:,:,:]= -deriv(P[:,:,:],i)

figure(6)
pcolor(xx,yy,P[1,:,:])
xlabel('x')
ylabel('y')
axis("equal")

#Include here the PDF of the pressure

########################################################
##ASSIGNMENT 7: Calculate the rate-of-strain tensor    ##
########################################################

STR=zeros((3,3,N,N,N))

for i in range(0,3):
    for j in range(0,3):
        STR[i,j,:,:,:]=0.5*(GU[i,j,:,:,:] + GU[j,i,:,:,:])
        
########################################################
##ASSIGNMENT 8: Calculate the rate-of-rotation tensor  ##
########################################################

OME=zeros((3,3,N,N,N))

for i in range(0,3):
    for j in range(0,3):
        OME[i,j,:,:,:]=0.5*(GU[i,j,:,:,:] - GU[j,i,:,:,:])

        
########################################################
##ASSIGNMENT 9: Calculate the visc. forces and its div.## 
########################################################

#Viscous forces (vector)
VF=zeros((3,N,N,N))

#Loop over three directions
for i in range(0,3):
    for j in range(0,3):
        VF[i,:,:,:]=VF[i,:,:,:] + 2.0*nu*deriv(STR[i,j,:,:,:],j)  
        
        
#The divergence of the viscous forces (scalar)
DVF=zeros((N,N,N))        

#Loop over three directions
for i in range(0,3):
    DVF[:,:,:]=DVF[:,:,:] + deriv(VF[i,:,:,:],i)
     

########################################################
##ASSIGNMENT 9: Check that NL are conservative         ## 
########################################################
#Multiply the NL forces by u and take the mean

#NL is a vector with the non-linear in Navier--Stokes
#U is a vector
#The product of U with NL is ENL, which is a scalar (inner product)
ENL=zeros((N,N,N))

#Loop over three directions
for i in range(0,3):
        ENL[:,:,:]=ENL[:,:,:] + U[i,:,:,:]*(-NL[i,:,:,:] + GP[i,:,:,:])
    
#Take the mean of ENL and check that it is zero: Non-linear terms are conservative    
mean(ENL)
        
        
#########################################################
##ASSIGNMENT 9: Calculate the mean energy dissipation  ## 
#########################################################
#Multiply first the viscous forces by u and take the mean
#Use the formula -2nuS_ijS_ij and check that the mean is
#same

#First definition
DIS1=zeros((N,N,N))

#Second definition
DIS2=zeros((N,N,N))

for i in range(0,3):
        DIS1[:,:,:]=DIS1[:,:,:] + U[i,:,:,:]*(VF[i,:,:,:])
    
for i in range(0,3):
    for j in range(0,3):
        DIS2[:,:,:]=DIS2[:,:,:] - 2.0*nu*STR[i,j,:,:,:]*STR[i,j,:,:,:]

#Calculate the mean dissipation
epsilon=mean(DIS2)

########################################################
##ASSIGNMENT 10: Cal. viscous scale, time and velocity ## 
########################################################
##Calculate the characteristic variables of the small scales
##and compare them with the large scale. Get the Reynolds number
#and check that epsilons\sim urms^3/L

#Viscous scales

#Length scale
eta=(nu**3.0/abs(epsilon))**(1.0/4.0)

#Time scale
tau=(eta**2.0/abs(epsilon))**(1.0/3.0)

#Velocity scale
uta=eta/tau

#Large scales

#Size of the computational box
L=2.0*pi

#Velocity fluctuations
urms=sqrt(1.0/3.0*mean(U[0,:,:,:]**2.0 + U[1,:,:,:]**2.0 + U[2,:,:,:]**2.0))

#Time scale
T=L/urms

#Calculate two different Reynolds numbers
#Based on urms, L  and nu
Re_L=L*urms/nu

#Based on the time scales of large and small scales
Re_tau=(T/tau)**2.0

