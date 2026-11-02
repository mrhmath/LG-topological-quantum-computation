#%% Preamble
from IPython.display import display, Math, Latex
from sympy import * #load symbolic library
from sympy.polys.matrices import DomainMatrix
from math import comb
from itertools import product

q  = Symbol('q')
t  = Symbol('t') #t=q^beta
s  = Symbol('s') #s=q^alpha
a  = Symbol('alpha')
b  = Symbol('beta')

# choose maximum n to consider (cable index/tensorfactors-1)
bign = 2

def sparseM(x):
    return  ImmutableSparseMatrix(x)

def KP(*factors):
    prod = Matrix([[1]])
    for x in factors:
        prod = kronecker_product(prod,x)
    return  sparseM(prod)

def DMsm(x,a):
    return (x)*((a)*eye(shape(x)[1])).to_DM()

def qnum(k):
    return (q**k-q**(-k))*(q-q**(-1))**(-1)

def qfact(k):
    if k == 1 or k == 0:
        return 1
    return qnum(k)*qfact(k-1)

def qexp(x,k):
    if k==0:
        return x**0
    return qexp(x,k-1)+x**k*q**(comb(k,2))/qfact(k)

def swapst(x):
    return x.subs({s: t, t: s}, simultaneous=True)

def subsst(x):
    return simplify(x.subs([(a,log(s)/log(q)),(b,log(t)/log(q))]))

def sparsemult(a,b):
    return (a.to_DM()*b.to_DM()).to_Matrix()

#%% Matrix Definitions
# each generator X is specified on the representation V(n,alpha) by "X(n)" or "X(n,a)"

# setup algorithm for determining matrix entries
def F1gen(n,i,j):
    if j%4 == 0 and i==j+4 and (j//4)+1<=n:
        return 1
    if j%4 == 1 and i==j+1:
        return 1
    if j%4 == 2 and i==j+4 and (j//4)+1<=n:
        return qnum(2)
    if j%4 == 2 and i==j+2*4-1 and (j//4)+2<=n:
        return -1
    if j%4 == 3 and i==j+4 and (j//4)+1<=n:
        return 1
    else:
        return 0

def F2gen(n,i,j):
    if j%4 == 0 and i==j+1:
        return 1
    if j%4 == 2 and i==j+1:
        return 1
    else:
        return 0

def E1gen(n,i,j):
    if j%4 == 0 and i==j-4:
        return qnum(j//4)*qnum(n-(j//4)+1)
    if j%4 == 1 and i==j-4:
        return qnum(j//4)*qnum(n-(j//4)+1)
    if j%4 == 2 and i==j-1:
        return qnum(n+1-2*(j//4))
    if j%4 == 2 and i==j-4 and (j//4)-1>=0:
        return qnum(j//4)*qnum(n-(j//4)+1)
    if j%4 == 3 and i==j-4 and (j//4)-1>=0:
        return qnum(j//4)*qnum(n-(j//4)+1)
    else:
        return 0

def E2gen(n,x,i,j):
    if j%4 == 1 and i==j-1:
        return qnum(x+(j//4))
    if j%4 == 2 and i==j+4-2 and (j//4)+1<=n:
        return qnum(x+(j//4))
    if j%4 == 3 and i==j-1:
        return qnum(x+(j//4)+1)
    if j%4 == 3 and i==j+4-2 and (j//4)+1<=n:
        return -qnum(x+(j//4))
    else:
        return 0

def K1gen(n,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return q**(n-2*(j//4))
    if j%4 == 1:
        return q**(n-2*(j//4)+1)
    if j%4 == 2:
        return q**(n-2*(j//4+1)+1)
    if j%4 == 3:
        return q**(n-2*(j//4+1)+2)

def K2gen(n,x,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return q**(x+(j//4))
    if j%4 == 1:
        return q**(x+(j//4))
    if j%4 == 2:
        return q**(x+((j//4)+1))
    if j%4 == 3:
        return q**(x+((j//4)+1))

def H1gen(n,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return (n-2*(j//4))
    if j%4 == 1:
        return (n-2*(j//4)+1)
    if j%4 == 2:
        return (n-2*(j//4+1)+1)
    if j%4 == 3:
        return (n-2*(j//4+1)+2)

def H2gen(n,x,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return (x+(j//4))
    if j%4 == 1:
        return (x+(j//4))
    if j%4 == 2:
        return (x+((j//4)+1))
    if j%4 == 3:
        return (x+((j//4)+1))

def bosgen(n,i,j):
    if i!=j:
        return 0
    if j%4 == 0:
        return (-1)**n
    if j%4 == 1:
        return (-1)**(n+1)
    if j%4 == 2:
        return (-1)**(n+1)
    if j%4 == 3:
        return (-1)**n

# function to output the matrix
def F1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: F1gen(n,i,j)))

def F2(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: F2gen(n,i,j)))

def F12(n):
    return F2(n)*F1(n)-q**-1*F1(n)*F2(n)

def E1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: E1gen(n,i,j)))

def E2(n,x):
    return simplify(ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: E2gen(n,x,i,j)))

def E12(n,x):
    return E1(n)*E2(n,x)-q*E2(n,x)*E1(n)

def K1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: K1gen(n,i,j)))

def K2(n,x):
    return simplify(ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: K2gen(n,x,i,j)))

def H1(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: H1gen(n,i,j)))

def H2(n,x):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: H2gen(n,x,i,j)))

def bos(n):
    return (ImmutableSparseMatrix(4*n+4, 4*n+4, lambda i, j: bosgen(n,i,j)))

#%% Coproducts

def DF1(n,m):
    return KP(F1(n), eye(4*m+4)) + KP(K1(n)**(-1), F1(m))

def DF2(n,m,x):
    return KP(F2(n), eye(4*m+4)) + KP(K2(n,x)**(-1)*bos(n), F2(m))

def DE1(n,m):
    return KP(E1(n), K1(m)) + KP(eye(4*n+4), E1(m))

def DE2(n,m,x,y):
    return KP(E2(n,x), K2(m,y)) + KP(bos(n), E2(m,y))

def DK1(n,m):
    return KP(K1(n),K1(m))

def DK2(n,m,x,y):
    return KP(K2(n,x),K2(m,y))

def DH1(n,m):
    return KP(H1(n), eye(4*m+4)) + KP(eye(4*n+4), H1(m))

def DH2(n,m,x,y):
    return KP(H2(n,x), eye(4*m+4)) + KP(eye(4*n+4), H2(m,y))

def Dbos(n,m):
    return KP(bos(n),bos(m))

# DM and t,s variable version
def DMDF1(n,m):
    return (DF1(n,m)).to_DM()

def DMDF2(n,m):
    return simplify(DF2(n,m,log(t)/log(q))).to_DM()

def DMDE1(n,m):
    return (DE1(n,m)).to_DM()

def DMDE2(n,m):
    return (simplify(DE2(n,m,log(t)/log(q),log(s)/log(q)))).to_DM()

def DMDK1(n,m):
    return (DK1(n,m)).to_DM()

def DMDK2(n,m):
    return simplify(DK2(n,m,log(t)/log(q),log(s)/log(q))).to_DM()

def DMDH1(n,m):
    return (DH1(n,m)).to_DM()

def DMDH2(n,m):
    return simplify(DH2(n,m,log(t)/log(q),log(s)/log(q))).to_DM()

def DMDbos(n,m):
    return Dbos(n,m).to_DM()

#%% R-matrix/braiding

def Pgen(n,m,i,j):
    if i//(4*(n+1)) == j%(4*(m+1)) and i%(4*(n+1)) == j//(4*(m+1)):
        return 1
    else:
        return 0

def P(n,m):
    return (Matrix(16*(n+1)*(m+1), 16*(n+1)*(m+1), lambda i, j: Pgen(n,m,i,j))).to_DM()

# multiplies the factor q^(2ab) to simplify the handling of the computation
def expH(n,m):
    H2n = H2(n,log(t)/log(q))
    H2m = H2(m,log(s)/log(q))
    exph1 = simplify(-(2*KP(H2n, H2m)+KP(H2n,H1(m))+KP(H1(n),H2m)).diagonal())
    exph = [q**(i+2*log(s)*log(t)/log(q)**2) for i in exph1]
    eH = zeros(16*(n+1)*(m+1))
    for i in range(16*(n+1)*(m+1)):
        eH[i,i]=exph[i]
    return ImmutableSparseMatrix(simplify(eH)).to_DM()


def quasiR(n,m):
    ID = KP(eye(4*(n+1)),eye(4*(m+1))).to_DM()
    DME1F1 = ((q-q**-1) * KP(E1(n),F1(m)))#.to_DM()
    exp1   = qexp(DME1F1 , min(n,m)+1).to_DM()
    exp12  = ID - ((q-q**-1) * KP(E12(n,log(t)/log(q))*bos(n) , F12(m))).to_DM()
    exp2   = ID - ((q-q**-1) * KP( E2(n,log(t)/log(q))*bos(n) , F2(m) )).to_DM()
    return (exp2*exp12*exp1)

def Bos(n,m):
    return ((KP(eye(4*(n+1)),eye(4*(m+1)))+KP(bos(n),eye(4*(m+1)))+KP(eye(4*(n+1)),bos(m))-KP(bos(n),bos(m)))/2).to_DM()

def DMR(n,m):
    return P(n,m)*Bos(n,m)*expH(n,m)*quasiR(n,m)

def R(n,m):
    return simplify(DMR(n,m).to_Matrix())


x=symbols('x0:4')
y=symbols('y0:4')
R00=R(0,0)

def totensor(x):
    return (x//4,x%4)

def condense(R):
    C = zeros(4)
    for j in range(16):
        for i in range(16):
            xind,yind = totensor(i)
            entry = factor(R[i,j])
            C[totensor(j)] += expand(numer(entry))/expand(denom(entry))*x[xind]*y[yind]
    return C

display(Math(r"q^{2\alpha\beta}(c_{V,W}(x_{i,j}))="),condense(R00))

#%% Highest Weight Vectors in V(0,b)xV(n,a)
# decomposes with highest weight vectors hstay(s), hright(u), hup(r), hleft(l)

def DE212F212(n):
    return (DMDE2(0,n)*DMDE1(0,n)*DMDE2(0,n)*DMDF2(0,n)*DMDF1(0,n)*DMDF2(0,n))


fullEF=[]
for i in range(bign+1):
    fullEF.append(DE212F212(i))

def hstay(n):
    return KP(eye(4)[:,0], eye(4*(n+1))[:,0]).to_DM()

def hright(n):
    if n<bign+1:
        w0 = fullEF[n]*KP(eye(4)[:,0], F2(n)[:,0]).to_DM()
        w0m = w0.to_Matrix()
        return DMsm(w0,w0m[next((index for index,value in enumerate(w0m) if value != 0), None)]**-1)
    else:
        return DE212F212(n,x,y)*KP(eye(4)[:,0], F2(n)[:,0]).to_DM()

def hup(n):
    if n<bign+1:
        w0 = fullEF[n]*KP(eye(4)[:,0], F2(n)*F1(n)*F2(n)[:,0]).to_DM()
        w0m = w0.to_Matrix()
        return DMsm(w0,w0m[next((index for index,value in enumerate(w0m) if value != 0), None)]**-1)
    else:
        return DE212F212(n,x,y)*KP(eye(4)[:,0], F2(n)*F1(n)*F2(n)[:,0]).to_DM()

def hleft(n):
    if n<bign+1:
        w0 = fullEF[n]*KP(eye(4)[:,0], qnum(n)*F1(n)*F2(n)[:,0]-qnum(n+1)*F2(n)*F1(n)[:,0]).to_DM()
        w0m = w0.to_Matrix()
        return DMsm(w0,w0m[next((index for index,value in enumerate(w0m) if value != 0), None)]**-1)
    else:
        return DE212F212(n)*KP(eye(4)[:,0],qnum(n)*F1(n)*F2(n)[:,0]-qnum(n+1)*F2(n)*F1(n)[:,0]).to_DM()


#%% Change of basis matrix for V(0,b)xV(n,a) tensor decomp

def basismat(n):
    ID = KP(eye(4),eye(4*(n+1))).to_DM()
    DF2Mat = DMDF2(0,n)
    DF1Mat = DMDF1(0,n)
    DE2Mat = DMDE2(0,n)
    DE1Mat = DMDE1(0,n)
    #generated by hstay V(n,a+b)
    w0 = hstay(n)
    B = w0
    for i in range(n+1):
        w = w0
        for J in [ID, DF2Mat, DF1Mat, DF2Mat]:
            w = (J*w)
            B = B.hstack(w)
        w0 = DF1Mat*w0
    # delete repeated first entry from the construction (cannot initialize DM)
    B = B[:,1:]
    #generated by hup V(n,a+b+1)
    w0 = hup(n)
    w0m = w0.to_Matrix()
    w0 = DMsm(w0,w0m[next((index for index,value in enumerate(w0m) if value != 0), None)]**-1)
    for i in range(n+1):
        w = w0
        for J in [ID, DF2Mat, DF1Mat, DF2Mat]:
                w = J*w
                B = B.hstack(w)
        w0 = DF1Mat*w0
    #generated by hright V(n+1,a+b)
    w0 = hright(n)
    w0m = w0.to_Matrix()
    w0 = DMsm(w0,w0m[next((index for index,value in enumerate(w0m) if value != 0), None)]**-1)
    for i in range(n+2):
        w = w0
        for J in [ID, DF2Mat, DF1Mat, DF2Mat]:
            w = J*w
            B = B.hstack(w)
        w0 = DF1Mat*w0
    #generated by hleft V(n-1,a+b+1)
    if n>0:
        w0 = hleft(n)
        w0m = w0.to_Matrix()
        w0 = DMsm(w0,w0m[next((index for index,value in enumerate(w0m) if value != 0), None)]**-1)
        for i in range(n):
            w = w0
            for J in [ID, DF2Mat, DF1Mat, DF2Mat]:
               w = J*w
               B = B.hstack(w)
            w0 = DF1Mat*w0
    return B

#%% Change of basis and forks

bases = []
basesinv = []
for i in range(bign):
    DMbi = (basismat(i))
    DMbiinv = DMbi.inv()
    bases.append(DMbi.to_Matrix())
    basesinv.append(DMbiinv.to_Matrix())

# V is upward fork: V('Vcolor','Vind') -> V(0,lcolor) x V(xind, rcolor)
# xdir = direction: stay, up, right, left
# xind = index n
def V(xdir,xind, lcolor,rcolor):
        if xdir == 's':
            return bases[xind][:,0:4*(xind+1)].subs({s: rcolor, t: lcolor}, simultaneous=True)
        if xdir == 'u':
            return bases[xind][:,4*(xind+1):8*(xind+1)].subs({s: rcolor, t: lcolor}, simultaneous=True)
        if xdir == 'r':
            return bases[xind][:,8*(xind+1):8*(xind+1)+4*(xind+2)].subs({s: rcolor, t: lcolor}, simultaneous=True)
        if xdir == 'l':
            return bases[xind][:,12*xind+16:12*xind+16+4*(xind)].subs({s: rcolor, t: lcolor}, simultaneous=True)

# L is downward fork: V(0,lcolor) x V(xind, rcolor) -> V(n_x,a_x)
def L(xdir,xind, lcolor,rcolor):
        if xdir == 's':
            return basesinv[xind][0:4*(xind+1),:].subs({s: rcolor, t: lcolor}, simultaneous=True)
        if xdir == 'u':
            return basesinv[xind][4*(xind+1):8*(xind+1),:].subs({s: rcolor, t: lcolor}, simultaneous=True)
        if xdir == 'r':
            return basesinv[xind][8*(xind+1):8*(xind+1)+4*(xind+2),:].subs({s: rcolor, t: lcolor}, simultaneous=True)
        if xdir == 'l':
            return basesinv[xind][12*xind+16:12*xind+16+4*(xind),:].subs({s: rcolor, t: lcolor}, simultaneous=True)


# Vcolor determines color through a fork
def Vcolor(xdir, xind, lcolor,rcolor):
    if xdir == 's' or xdir == 'r':
        return lcolor*rcolor
    if xdir == 'l' or xdir == 'u':
        return q*lcolor*rcolor

# Vind determines index through a fork
def Vind(xdir, xind, lcolor,rcolor):
    if xdir == 's' or xdir == 'u':
        return xind
    if xdir == 'r':
        return xind+1
    if xdir == 'l':
        return xind-1

# Computes the right branching fork V colored by s on the top right and t on the left
# the directions from the top are w = ['x1dir', 'x2dir', ...]
def bigV(w):
    stack = eye(4)
    stackcolor = s
    ind = 0
    for x in w:
        stack = sparsemult(KP(eye(4),stack),V(x,ind, t, stackcolor))
        stackcolor = Vcolor(x,ind, t, stackcolor)
        ind = Vind(x,ind, t, stackcolor)
    return stack

# Computes the right branching fork L colored by s on the bottom
# the directions from the bottom are w = ['x1dir', 'x2dir', ...]
def bigL(w):
    stack = eye(4)
    stackcolor = s
    ind = 0
    for x in w:
        stack = sparsemult(L(x,ind, t, stackcolor),KP(eye(4),stack))
        stackcolor = Vcolor(x,ind, t, stackcolor)
        ind = Vind(x,ind, t, stackcolor)
    return stack

names = ['su','us','rl']
name_to_index = {'su':0, 'us':1, 'rl':2}
name_to_path = {'su':['s','u'], 'us':['u','s'], 'rl':['r','l']}

name_to_V={}
name_to_L={}
for i in names:
    name_to_V[i] = bigV(name_to_path[i])
    name_to_L[i] = bigL(name_to_path[i])

R00eye = KP(R00.subs(s,t),eye(4))
eyeR00sq = KP(eye(4),swapst(R00)*R00)
def mat_2_components(mat):
    output = zeros(3,3)
    for i in names:
        for j in names:
            M = simplify((name_to_L[j]*mat*name_to_V[i])).applyfunc(factor)
            assert simplify(M-M[0,0]*eye(4))==zeros(4,4)
            output[name_to_index[j],name_to_index[i]] = M[0,0]
    return output


# divides out the earlier factor q^(2ab) with b=-1/2
Remat = mat_2_components(s*R00eye).subs(t,q**Rational(-1,2))
display(Math(r"\rho_1(\sigma_1)="),Remat)
quo, rem = div(expand(numer(Remat[2,2])),expand(denom(Remat[2,2])),s)
display(Math(r"\rho_1(\sigma_1)_{rl,rl}="),quo+factor(rem)/factor(expand(denom(Remat[2,2]))))
eR2mat = mat_2_components(s**2*eyeR00sq).subs(t,q**Rational(-1,2))
display(Math(r"\rho_1(J_2)="),eR2mat)

