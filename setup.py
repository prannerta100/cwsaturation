import os
import numpy as np
from scipy.sparse import csr_matrix, spdiags, coo_matrix
from scipy.ndimage import gaussian_filter1d
from scipy.sparse import block_diag, identity,bmat,diags, spdiags
from scipy.sparse.linalg import gmres, spsolve
from lmfit import minimize, Parameters, report_fit #for pp vs b1
#from scikits.umfpack import spsolve
#from pypardiso import spsolve
import time
#from arnoldi import *
from math import ceil
import copy
import numpy as np
import matplotlib.pyplot as plt

#Hyperparameters
MIN_LEMX = 8#16
MXDIM = 2000
MXEL = 300*MXDIM
KRYLOV_THRESH = 2000#0
MXPT = 4096

# modify limits.inc; going for a small MXDIM as X band non-ultra-slow motions
f = open('limits.inc','w')
f.write('''
    integer :: MXDIM,MXEL,MXLVAL,NAMELG,MMAX,MXPT
    parameter (MXDIM={},MXEL={},MXLVAL=200,NAMELG=40,MMAX=100,MXPT={})
'''.format(MXDIM, MXEL, MXPT))
f.close()

#install the "mat" package from the Fortran90 files listed below using the Numpy f2py wrapper; this is Windows version!!!!
#see "quick and smart way" at https://numpy.org/doc/stable/f2py/f2py.getting-started.html#the-quick-and-smart-way
#ignore the warnings generated during numpy.f2py; have a lot of unused variables; don't know why the "uninitialized-maybe" warnings come 
print('No need to compile the Fortran code again and again, go to setup.py and comment the f2py line to avoid wasting time every compilation')
#os.system('f2py -c --fcompiler=gnu95 --compiler=mingw32 -m mat  generate_matrices.f90 stveco.f90 anxlk.f90 matrxd.f90 matrxo.f90 cd2km.f90 fz.f90 ccrint_new.f90 bessel.f90 ipar.f90 plgndr.f90 w3j.f90')
os.system('python3 -m numpy.f2py -c -m mat generate_matrices.f90 stveco.f90 anxlk.f90 matrxd.f90 matrxo.f90 cd2km.f90 fz.f90 ccrint_new.f90 bessel.f90 ipar.f90 plgndr.f90 w3j.f90')
os.system('python3 -m numpy.f2py -c -m gconvl gconvl.f')

import mat #the Fortran90 wrapped module just generated, mat.generate_matrices is a F90 subroutine that can be called in Python 
import gconvl
# definitions of parameters as per their positions in the arrays fed to Fortran
params_double_def = {'gxx':0, 'gyy':1, 'gzz':2, 'axx':3, 'ayy':4, 'azz':5, 'dx':6, 'dy':7, 'dz':8,\
                     'pml':9, 'pmxy':10, 'pmzz':11, 'djf':12, 'djfprp':13, 'oss':14, 'psi':15,\
                     'ald':16, 'bed':17, 'gad':18, 'alm':19, 'bem':20, 'gam':21,\
                     'c20':22, 'c22':23, 'c40':24, 'c42':25, 'c44':26,\
                     't2edi':27, 't2ndi':28, 't2efi':29, 't1edi':30, 't1ndi':31,\
                     'b0':32, 'a0':33, 'g0':34, 'pl':35, 'pkxy':36, 'pkzz':37}
params_int_def = {'in2':0, 'ipdf':1, 'ist':2,\
                  'lemx':3, 'lomx':4, 'kmx':5, 'mmx':6, 'ipnmx':7}

# suggest a decent basis size as per Budil's equation; for fast motions best is to use lemx=8 at X band
def budil_basis_size(rbar, B0):
    '''
    Input:  rbar-> eztimate of the motional rate
            B0 -> static magnetic field (Gauss)
    Output: tuple of truncation indices (lemx,lomx,kmx,mmx) as per Budil's empirical equation presented in an ACERT workshop
    '''
    lemx=round((3.07+0.0068*0.0028*B0)*(8.40+0.0075*0.0028*B0-rbar+0.3010)**2)
    lemx=MIN_LEMX if lemx<MIN_LEMX else lemx
    lemx+=lemx%2
    lomx=round(0.7*lemx)
    lomx+=lomx%2+1
    kmx=round(lemx*0.5)
    kmx+=kmx%2
    mmx=kmx
    return lemx,lomx,kmx,mmx

def conv2coo(filename):
    '''
    misc function that converts <filename> to a scipy.sparse.COO matrix; 
    assuming the file comes from a language like Fortran/Matlab that begins indexing from 1
    to convert to Python we need to subtract indices by 1, b/c Python starts arrays from 0

    Input: filename string
    Output: COO matrix
    '''
    print('Reading ' + filename)
    x = np.loadtxt(filename)
    dim = int(max(x[:,0]))
    I = x[:,0].astype(int) - 1
    J = x[:,1].astype(int) - 1
    print('Dimensions '+str(dim)+'x'+str(dim))
    E = x[:,2]
    return coo_matrix((E,(I,J)),shape=(dim,dim))

def generate_from_params(offdiag_basis_file, simparams_double, simparams_int):
    '''
    Input:
    1. offdiag_basis_file: input basis set (a nonsensical, non-existent file like 'xoxo' means the F90 subroutine 
    will assume that you are generating matrices, etc. from scratch) 
    2. simparams_double: simulation parameters (look Budil et al 1996 Table 1 to find the meanings of each parameter
    gxx(1); gyy(2); gzz(3);axx(4); ayy(5); azz(6);dx(7); dy(8); dz(9);pml(10); pmxy(11); pmzz(12);djf(13); djfprp(14);oss(15); psi(16);ald(17); bed(18); gad(19);alm(20); bem(21); gam(22);
    c20(23); c22(24);c40(25); c42(26); c44(27);t2edi(28);t2ndi(29); t2efi(30);t1edi(31); t1ndi(32);b0(33); a0(34); g0(35);pl(36); pkxy(37);pkzz(38);
    3. simparams_int: more simulation parameters that are integers, not double precision numbers 
    in2(1); ipdf(2); ist(3); lemx(4);lomx(5); kmx(6); mmx(7); ipnmx(8)

    Output: 
    matx: off-diag space matrix in scipy.csr format
    matz: diag space matrix in scipy.csr format
    pp: pulse propagator in scipy.csr format
    stvx: off-diag space starting vector
    
    Notes:
    1. <ndimo_in> has been removed, PG 2/17/21 
    2. To people new to SLE, off-diagonal and diagonal spaces are DIFFERENT from off-diagonal and diagonal parts of a matrix. 
    At a deeper level there is a relation, but you need to know density matrices in quantum mechanics.
    '''
    #print(simparams_int)
    zmat_offdiag, zdiag_offdiag, izmat_offdiag, jzmat_offdiag, kzmat_offdiag, zmat_diag, zdiag_diag, \
    izmat_diag, jzmat_diag, kzmat_diag, mpid, mpp, stvx, nelreo, nelimo, ndimo, \
    nelred, nelimd, ndimd = mat.generate_matrices(offdiag_basis_file,simparams_double, simparams_int) #'xoxo' is a dummy file
    #print(nelreo, nelimo, ndimo, nelred, nelimd, ndimd) 
    #off-diag space starting vector <stvx>
    #such truncations viz. [:ndimo] are needed, because F90 f2py doesn't allow allocatable arrays it seems, so we need to define a much larger array
    #https://stackoverflow.com/questions/34579769/f2py-error-with-allocatable-arrays/34708146
    stvx=stvx[:ndimo] 
    #off-diag space SLE matrix <matx>
    offi=csr_matrix((zmat_offdiag[:nelimo],izmat_offdiag[:nelimo]-1, jzmat_offdiag[:(ndimo+1)]-1),shape=(ndimo,ndimo))
    if nelreo > 0:
        offr=csr_matrix((zmat_offdiag[:(-nelreo-1):-1],izmat_offdiag[:(-nelreo-1):-1]-1,kzmat_offdiag[:(ndimo+1)]-1),shape=(ndimo,ndimo))
    else:
        offr=csr_matrix((ndimo,ndimo))
    matx=offr+offr.transpose()+spdiags(zdiag_offdiag[0,:ndimo],[0],ndimo,ndimo)+1.0j*(offi+offi.transpose()+spdiags(zdiag_offdiag[1,:ndimo],[0],ndimo,ndimo))
    #diag space SLE matrix <matz>
    offi=csr_matrix((zmat_diag[:nelimd],izmat_diag[:nelimd]-1, jzmat_diag[:(ndimd+1)]-1),shape=(ndimd,ndimd))
    if nelred > 0:
        offr=csr_matrix((zmat_diag[:(-nelred-1):-1],izmat_diag[:(-nelred-1):-1]-1,kzmat_diag[:(ndimd+1)]-1), shape=(ndimd,ndimd))
    else:
        offr=csr_matrix((ndimd,ndimd))
    #print(offr.shape, spdiags(zdiag_diag[0,:ndimd],[0],ndimd,ndimd).shape, spdiags(zdiag_diag[1,:ndimd],[0],ndimd,ndimd).shape, offi.shape)
    matz=offr+offr.transpose()+spdiags(zdiag_diag[0,:ndimd],[0],ndimd,ndimd)+1.0j*(offi+offi.transpose()+spdiags(zdiag_diag[1,:ndimd],[0],ndimd,ndimd))
    #pulse propagator <pp>
    mpp=mpp[:ndimd]
    mpid=mpid[:ndimo]
    indx=[]
    for k in range(ndimo):
        if mpid[k]==1:
            indx.append(k)
        elif mpid[k]==2:
            indx.append(k)
            indx.append(k)
    pp=coo_matrix((mpp,(np.arange(ndimd),np.array(indx))),shape=(ndimd,ndimo)).tocsr()/np.sqrt(2)
    return matx, matz, pp, stvx #, [zmat_offdiag, zdiag_offdiag, izmat_offdiag, jzmat_offdiag, kzmat_offdiag],\
                                # [zmat_diag, zdiag_diag, izmat_diag, jzmat_diag, kzmat_diag]


def cw_spec(bgrid=np.linspace(-60,60,128)+3360, params_in=dict(), basis_file='xoxo', prune_on=0):
    '''
    calculates the derivative spectrum for a given magnetic field grid, basis file input
    
    Inputs
    ------
    bgrid: grid of magnetic field values in Gauss, need not be uniformly spaced
    params_in: dictionary of parameters
    basis_file: input basis file; very unlikely this will be used for saturation calculations
    prune_on: integer; 0 means no prune, 1 means prune matx, use the pruned matx to prune matz and then proceed
    
    Output
    ------
    tuple of bgrid  and the derivative spectrum calculated by forward difference
    bgrid is an input parameter, so redundant to output bgrid; will change in a future version    
    '''
    simparams_double=np.array(([2.008820, 2.006200, 2.002330, 5.20,   5.80,  34.40, 8.18, 8.18, 9.27, 0,0,0,0,0, 0.0, 45, 0,0,0,0,0,0,2.0, 0.0, 0,0,0,0.0,0,0,np.log10(2*8.8e4),0,3360,0,0,0,0,0]))
    lemx, lomx, kmx, mmx = budil_basis_size(simparams_double[params_double_def['dx']],simparams_double[params_double_def['b0']]) #[12,9,4,4]
    simparams_int=np.array([2,0,0,lemx,lomx,kmx,mmx,2])#([2,0,0,22,13,7,7,2])#([2,0,0,44,33,14,14,2])
    #simparams_int=np.array([2,0,0,22,19,14,2,2])
    # read parameters from the dictionary
    for x in params_in:
        if x in params_double_def:
            simparams_double[params_double_def[x]] = params_in[x]
        if x in params_int_def:
            simparams_int[params_int_def[x]] = params_in[x]
    # off-diagonal space shift shiftx (same as lb!)
    shiftx = params_in['shiftx'] if 'shiftx' in params_in else 0.0
    # diagonal space shift shiftz
    shiftz = params_in['shiftz'] if 'shiftz' in params_in else 0.0
    # prune tol
    ptol = params_in['ptol'] if 'ptol' in params_in else 0.0 #001
    # gmres tol
    gmres_tol = params_in['gmres_tol'] if 'gmres_tol' in params_in else 0.0000001
    # overall scaling factor
    scale = params_in['scale'] if 'scale' in params_in else 1.0
    # overall x axis shift factor
    shiftg = params_in['shiftg'] if 'shiftg' in params_in else 0.0
    # gib0
    gib0 = params_in['gib0'] if 'gib0' in params_in else 0.0
    # gib2
    gib2 = params_in['gib2'] if 'gib2' in params_in else 0.0
    # nort 
    nort = int(params_in['nort']) if 'nort' in params_in else 10
    # print parameters
    print(dict(zip(params_double_def.keys(),simparams_double)))
    print(dict(zip(params_int_def.keys(),simparams_int)))
    # b0, should be there in simparams_double
    B0 = simparams_double[params_double_def['b0']]
    # b1, should be there in params_in
    B1 = params_in['b1']
    print('Computing '+str(B1)+' Gauss')
    #cfact=1e-06*np.mean(simparams_double[:3])*9.2731e-21 / 1.05443e-27
    #omarrG=bgrid*2*np.pi/cfact
    omarrG=bgrid+shiftg-B0
    basis_file_trunc = 'xoxo'
    res=np.zeros_like(omarrG)
    #print('Computing '+str(B1)+' Gauss')
    # prune the off-diag space matrices; prune=1 means prune matx, use it to prune everything else
    # will add prune=2 for the case of pruning post mat_full creation
    if prune_on==1:
        ommin,ommax = -25,25
        prune_bgrid = np.linspace(ommin,ommax,20)
        simparams_double1=copy.deepcopy(simparams_double) #np.array([2.0084,2.0054,2.0019,5.0,5.0,32.6,5.3622,5.3622,6.6544,0,0,0,0,0,5.646,45,0,0,0,0,0,0,2.2572,-2.1782,0,0,0,6.733,0,0,5.568,0,6167.6,0,0,0,0,0])
        simparams_int1=copy.deepcopy(simparams_int) #np.array([2,0,0,lemx,lomx,kmx,mmx,2])#([2,0,0,22,13,7,7,2])#([2,0,0,44,33,14,14,2])
        simparams_double1[params_double_def['psi']] = 0.00001 #prune for one orientation
        matx1, matz1, pp1, stvx1 = generate_from_params(basis_file_trunc, simparams_double1, simparams_int1)
        matx1 += 1.0j*B0*identity(matx1.shape[0])
        prune_resv = np.zeros((matx1.shape[0],len(prune_bgrid)))
        for i in range(len(prune_bgrid)):
            m = matx1+(shiftx-1.0j*prune_bgrid[i]+1.0j*B0)*identity(matx1.shape[0])
            InvPrec = spdiags(1/m.diagonal(),[0],m.shape[0],m.shape[1])
            invec = spsolve(m,stvx1)
            prune_resv[:,i]=np.abs(invec/(stvx1.conjugate().transpose() @ invec))
        prune_offdiag = np.max(prune_resv,axis=1) > ptol
        prune_diag = (pp1 @ prune_offdiag) != 0
        # prune the offdiag matrix
        matx1=(matx1[prune_offdiag,:].tocsc())[:,prune_offdiag].tocsr()
        # prune the off-diag space starting vector
        stvx1=stvx1[prune_offdiag]
        # prune the diag space matrix
        matz1=(matz1[prune_diag,:].tocsc())[:,prune_diag].tocsr()
        # prune the pulse propagator
        pp1=(pp1[prune_diag,:].tocsc())[:,prune_offdiag].tocsr()

    if nort > 0: #MOMD
        for iort in range(nort):
            cspsi=iort/(nort-1) #epsilon to avoid psi=0 exactly
            gib = gib0 + gib2*(1-cspsi**2)
            wline = np.sqrt(gib*gib+shiftx*shiftx)
            if cspsi == 1:
                cspsi -= 1.0e-6
            simparams_double1=copy.deepcopy(simparams_double)#np.array([2.0084,2.0054,2.0019,5.0,5.0,32.6,5.3622,5.3622,6.6544,0,0,0,0,0,5.646,45,0,0,0,0,0,0,2.2572,-2.1782,0,0,0,6.733,0,0,5.568,0,6167.6,0,0,0,0,0])
            simparams_int1=copy.deepcopy(simparams_int)#np.array([2,0,0,lemx,lomx,kmx,mmx,2])#([2,0,0,22,13,7,7,2])#([2,0,0,44,33,14,14,2])
            simparams_double1[params_double_def['psi']] = np.arccos(cspsi)*180.0/np.pi
            print([simparams_double1])
            #print(simparams_int1)
            scal_momd = 0.5/(nort-1) if iort == 0 or iort == nort-1 else 1.0/(nort-1)
            matx1, matz1, pp1, stvx1 = generate_from_params(basis_file_trunc, simparams_double1, simparams_int1)
            matx1 += 1.0j*B0*identity(matx1.shape[0])
            if prune_on==1: # prune
                matx1=(matx1[prune_offdiag,:].tocsc())[:,prune_offdiag].tocsr()
                stvx1=stvx1[prune_offdiag]
                matz1=(matz1[prune_diag,:].tocsc())[:,prune_diag].tocsr()
                pp1=(pp1[prune_diag,:].tocsc())[:,prune_offdiag].tocsr()
            mat_full=bmat([[matx1, 0.5j*B1*pp1.transpose(), None],[0.5j*B1*pp1, matz1, -0.5j*B1*pp1],[None, -0.5j*B1*pp1.transpose(), matx1.conjugate().transpose()]])
            ndimo = matx1.shape[0]; ndimd=matz1.shape[0]
            stvx_full=np.hstack((1.0j*stvx1,np.zeros(ndimd),-1.0j*stvx1))
            stvx_full_left=abs(B1)*np.hstack((stvx1,np.zeros(ndimo+ndimd)))
            shifts=block_diag((shiftx*identity(ndimo),shiftz*identity(ndimd),shiftx*identity(ndimo)))
            signs=block_diag((identity(ndimo),0*identity(ndimd),-identity(ndimo)))
            '''
            mat_full=matx1
            stvx_full=stvx1
            stvx_full_left=abs(B1)*stvx1
            shifts=shiftx*identity(ndimo)
            signs=identity(ndimo)
            print(ndimo)
            '''
            tmpres = 0 * res
            if mat_full.shape[0] > KRYLOV_THRESH:
                for i in range(len(omarrG)):
                    InvPrec = diags(1/(mat_full+shifts-1.0j*omarrG[i]*signs).diagonal())
                    sol,info = gmres(mat_full+shifts-1.0j*omarrG[i]*signs,stvx_full,None,gmres_tol,200,ceil(mat_full.shape[0]/2000),InvPrec)
                    #sol,info = gmres(mat_full+shifts-1.0j*omarrG[i]*signs,stvx_full,None,gmres_tol,20,100,InvPrec)
                    if info > 0:
                        print("GMRES didn't converge for field offset "+str(omarrG[i])+", might be ok for other field values")
                    tmpres[i]=scal_momd*np.imag(stvx_full_left.transpose()@sol)                
            else:
                for i in range(len(omarrG)):
                    #sol = spsolve(mat_full+(shiftx-1.0j*omarrG[i])*identity(mat_full.shape[0]),stvx_full)
                    sol = spsolve(mat_full+shifts-1.0j*omarrG[i]*signs,stvx_full)
                    tmpres[i]=scal_momd*np.imag(stvx_full_left.transpose().conjugate()@sol)                #Transpose or Herm conj???????
            if wline > 0:            
                dummy_omarrG = np.linspace(min(omarrG), max(omarrG), 1000)
                #dummy_spec = np.sqrt(2*np.pi)*0.5*gaussian_filter1d(np.interp(dummy_omarrG, omarrG, tmpres), sigma=int(2*len(dummy_omarrG)*wline/(max(dummy_omarrG)-min(dummy_omarrG))))
                dummy_spec = gconvl.gconvl(np.hstack((np.interp(dummy_omarrG, omarrG, tmpres),np.zeros(MXPT-len(dummy_omarrG)))), wline, np.diff(dummy_omarrG)[0], 1000, 2048)[:len(dummy_omarrG)]
                res += np.interp(omarrG, dummy_omarrG, dummy_spec)
            else:
                res += tmpres
            '''
            for i in range(len(omarrG)):
                X = Q[:,:-1].transpose().conjugate() @ ((mat_full+shifts-1.0j*omarrG[i]*signs) @ Q[:,:-1])
                sol = Q[:,:-1] @ np.linalg.solve(X,np.eye(h.shape[1],1))
                res[i]+=np.real(1.0j*stvx_full_left.transpose().conjugate()@sol)
            '''
    else: #no MOMD
        print('nort was set to 0, will zero out psi and potential terms as well, no gib2 either')                    
        wline = np.sqrt(gib0*gib0+shiftx*shiftx)
        simparams_double1=copy.deepcopy(simparams_double)#np.array([2.0084,2.0054,2.0019,5.0,5.0,32.6,5.3622,5.3622,6.6544,0,0,0,0,0,5.646,45,0,0,0,0,0,0,2.2572,-2.1782,0,0,0,6.733,0,0,5.568,0,6167.6,0,0,0,0,0])
        simparams_int1=copy.deepcopy(simparams_int)#np.array([2,0,0,lemx,lomx,kmx,mmx,2])#([2,0,0,22,13,7,7,2])#([2,0,0,44,33,14,14,2])
        for x in ['c20','c22','psi']:
            simparams_double1[params_double_def[x]] = 0.0
        print([simparams_double1])
        matx1, matz1, pp1, stvx1 = generate_from_params(basis_file_trunc, simparams_double1, simparams_int1)
        matx1 += 1.0j*B0*identity(matx1.shape[0])
        if prune_on==1: # prune
            matx1=(matx1[prune_offdiag,:].tocsc())[:,prune_offdiag].tocsr()
            stvx1=stvx1[prune_offdiag]
            matz1=(matz1[prune_diag,:].tocsc())[:,prune_diag].tocsr()
            pp1=(pp1[prune_diag,:].tocsc())[:,prune_offdiag].tocsr()
        mat_full=bmat([[matx1, 0.5j*B1*pp1.transpose(), None],[0.5j*B1*pp1, matz1, -0.5j*B1*pp1],[None, -0.5j*B1*pp1.transpose(), matx1.conjugate().transpose()]])
        ndimo = matx1.shape[0]; ndimd=matz1.shape[0]
        stvx_full=np.hstack((1.0j*stvx1,np.zeros(ndimd),-1.0j*stvx1))
        stvx_full_left=abs(B1)*np.hstack((stvx1,np.zeros(ndimo+ndimd)))
        shifts=block_diag((shiftx*identity(ndimo),shiftz*identity(ndimd),shiftx*identity(ndimo)))
        signs=block_diag((identity(ndimo),0*identity(ndimd),-identity(ndimo)))
        tmpres = np.zeros_like(omarrG)
        for i in range(len(omarrG)):
            sol = spsolve(mat_full+shifts-1.0j*omarrG[i]*signs,stvx_full)
            tmpres[i] = np.imag(stvx_full_left.transpose()@sol)                
        # add wline        
        if wline > 0:
            dummy_omarrG = np.linspace(min(omarrG), max(omarrG), 1000)
            dummy_spec = gconvl.gconvl(np.hstack((np.interp(dummy_omarrG, omarrG, tmpres),np.zeros(MXPT-len(dummy_omarrG)))), wline, np.diff(dummy_omarrG)[0], 1000, 2048)[:len(dummy_omarrG)]
            #dummy_spec = np.sqrt(2*np.pi)*0.5*gaussian_filter1d(np.interp(dummy_omarrG, omarrG, tmpres), sigma=int(2*len(dummy_omarrG)*wline/(max(dummy_omarrG)-min(dummy_omarrG))))
            res = np.interp(omarrG, dummy_omarrG, dummy_spec)
        else:
            res = tmpres
    # return the derivative spectrum
    return bgrid,scale*np.gradient(res, omarrG)#np.hstack((0,np.diff(res)/np.diff(omarrG)))

