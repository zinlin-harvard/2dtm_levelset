#lenghts and distances are in units of wavelength
from dolfin import *

#domain and boundary labels
domext=0
domint=1
pmlx=2
pmly=3
pmlxy=5
bndDir=0
bndfacet=1
bndperx=2
bndpery=3
bndmat=4
bndNeu=5

#size of the entire domain including PMLs
Lx=1.0
Ly=1.0
Lpmlx=0.3
Lpmly=0.2
xleft=-0.3
ybottom=-0.2
xright=Lx+Lpmlx
ytop=Ly+Lpmly

#define wavelengths, frequencies and material parameters
lbda=0.4
freq=1/lbda
omega=2*pi*freq
omegasq=pow(omega,2)
eps0=complex(1.0,0.0)
eps1=complex(4.0,0.0)
eps0r=eps0.real
eps0i=eps0.imag
eps1r=eps1.real
eps1i=eps1.imag
mu0=complex(1.0,0.0)
mu1=complex(1.0,0.0)
muinv0r=(1.0/mu0).real
muinv0i=(1.0/mu0).imag
muinv1r=(1.0/mu1).real
muinv1i=(1.0/mu1).imag

#define the gaussian current position and width
Jcx=0.5
Jcy=0.5
Jrad=0.5

##############################################################################################
#mesh is imported from .xml file with pre-labelled subdomains
mesh=Mesh("mesh.xml")
domains=MeshFunction("size_t",mesh,"subdomains.xml")
boundaries=MeshFunction("size_t",mesh,"bnds.xml")

###############################################################################################

#Choose the finite element function space, here Lagrange elements of
#order 1 for scalar Hz field, for vector fields
#one should choose appropriate Nedelec curl elements 
V=FunctionSpace(mesh,"Lagrange",1)
W=V * V

#Universal simple Dirichlet boundary for all-around PML
def u0_boundary(x, on_boundary):
    return on_boundary
u0=(0.0,0.0)
bc = DirichletBC(W, u0, u0_boundary)

##########################################################################################
#Define expressions for anisotropic mu inverse pmls, ASSUMING vacuum index within pml
class MUINVPMLr(Expression):
    def __init__(self,xleft,xright,ytop,ybottom,Lpmlx,Lpmly,sigma,m):
        self.xleft=xleft
        self.xright=xright
        self.ytop=ytop
        self.ybottom=ybottom
        self.Lpmlx=Lpmlx
        self.Lpmly=Lpmly
        self.sigma=sigma
        self.m=m
    def eval(self,value,x):
        dx=min(x[0]-self.xleft,self.xright-x[0])
        dy=min(x[1]-self.ybottom,self.ytop-x[1])
        rhox=min(dx/self.Lpmlx,1.0)
        rhoy=min(dy/self.Lpmly,1.0)
        sximag=self.sigma*pow(1.0-rhox,self.m)
        syimag=self.sigma*pow(1.0-rhoy,self.m)
        valx=(1-sximag*sximag)/(pow(1-sximag*sximag,2)+4.0*pow(sximag,2))
        valy=(1-syimag*syimag)/(pow(1-syimag*syimag,2)+4.0*pow(syimag,2))
        value[0]=valx
        value[1]=0.0
        value[2]=0.0
        value[3]=valy
    def value_shape(self):
        return(2,2)

class MUINVPMLi(Expression):
    def __init__(self,xleft,xright,ytop,ybottom,Lpmlx,Lpmly,sigma,m):
        self.xleft=xleft
        self.xright=xright
        self.ytop=ytop
        self.ybottom=ybottom
        self.Lpmlx=Lpmlx
        self.Lpmly=Lpmly
        self.sigma=sigma
        self.m=m
    def eval(self,value,x):
        dx=min(x[0]-self.xleft,self.xright-x[0])
        dy=min(x[1]-self.ybottom,self.ytop-x[1])
        rhox=min(dx/self.Lpmlx,1.0)
        rhoy=min(dy/self.Lpmly,1.0)
        sximag=self.sigma*pow(1.0-rhox,self.m)
        syimag=self.sigma*pow(1.0-rhoy,self.m)
        valx=-2.0*sximag/(pow(1-sximag*sximag,2)+4.0*pow(sximag,2))
        valy=-2.0*syimag/(pow(1-syimag*syimag,2)+4.0*pow(syimag,2))
        value[0]=valx
        value[1]=0.0
        value[2]=0.0
        value[3]=valy
    def value_shape(self):
        return(2,2)

#pml parameters and epsinvpml tensor expressions
refl=1e-10
m=2
sigma=-(m+1)*ln(refl)/(2.0*omega*Lpmlx)
muinvpmlreal=MUINVPMLr(xleft,xright,ytop,ybottom,Lpmlx,Lpmly,sigma,m)
muinvpmlimag=MUINVPMLi(xleft,xright,ytop,ybottom,Lpmlx,Lpmly,sigma,m)

###########################################################################################################

# Redefine new measures for the domains
dx = Measure("dx")[domains]
#ds = Measure("ds")[boundaries]
dS = Measure("dS")[boundaries]

(ur,ui)=TrialFunctions(W)
(vr,vi)=TestFunctions(W)

#Electric current J position, radius and polarization
Jr=Expression('1/(pi*rad*rad)*exp(-(pow(x[0]-cx,2) + pow(x[1]-cy,2))/pow(rad,2))',cx=0.0,cy=0.0,rad=0.05)
Jr.cx=Jcx
Jr.cy=Jcy
Jr.rad=Jrad
Ji=Expression('value',value=0.0)

#define the common left hand side form
pmlgradreal=inner(nabla_grad(vr),muinvpmlreal*nabla_grad(ur)-muinvpmlimag*nabla_grad(ui))
pmlgradimag=inner(nabla_grad(vi),muinvpmlreal*nabla_grad(ui)+muinvpmlimag*nabla_grad(ur))
pmlgradform=pmlgradreal+pmlgradimag

mat0gradreal=inner(nabla_grad(vr),muinv0r*nabla_grad(ur)-muinv0i*nabla_grad(ui))
mat0gradimag=inner(nabla_grad(vi),muinv0r*nabla_grad(ui)+muinv0i*nabla_grad(ur))
mat0gradform=mat0gradreal+mat0gradimag

mat1gradreal=inner(nabla_grad(vr),muinv1r*nabla_grad(ur)-muinv1i*nabla_grad(ui))
mat1gradimag=inner(nabla_grad(vi),muinv1r*nabla_grad(ui)+muinv1i*nabla_grad(ur))
mat1gradform=mat1gradreal+mat1gradimag

epsform0real=omegasq*(vr*(eps0r*ur-eps0i*ui))
epsform0imag=omegasq*(vi*(eps0r*ui+eps0i*ur))
epsform0=epsform0real+epsform0imag
epsform1real=omegasq*(vr*(eps1r*ur-eps1i*ui))
epsform1imag=omegasq*(vi*(eps1r*ui+eps1i*ur))
epsform1=epsform1real+epsform1imag

pmlform=(-pmlgradform+epsform0)*(dx(pmlx)+dx(pmly)+dx(pmlxy))
mat0form=(-mat0gradform+epsform0)*dx(domext)
mat1form=(-mat1gradform+epsform1)*dx(domint)

a=pmlform + mat0form + mat1form

#define the right hand side form from J (electric current)
Jform= omega * (-vr * Ji + vi * Jr) * (dx(domint) + dx(domext) + dx(pmlx) + dx(pmly) + dx(pmlxy)) 
L=Jform

#assemble lhs and rhs and solve Aw = b
A=assemble(a)
solver=LUSolver(A)

b=assemble(L)
bc.apply(A,b)
w=Function(W)
solver.solve(w.vector(),b)
(Ezr,Ezi)=w.split()

import numpy as np

Efielddata=[]
grid = np.loadtxt("grid.txt")
for i in range(len(grid)):
    x=grid[i][0]
    y=grid[i][1]
    Efielddata.append([x,y,Ezr(Point(x,y)),Ezi(Point(x,y))])

np.savetxt("EzfieldGrid.txt",Efielddata)


