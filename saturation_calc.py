# -*- coding: utf-8 -*-
from setup import *
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
import gconvl

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


# boris no potential
'''
stt=time.time()
def sat_residual(params_fit=dict(), params_nonfit=dict(), bgrid=np.reshape(np.linspace(-60,60,256)+3360,(1,-1)), spec_expt=np.zeros((1,128)), b1_list = [0.1], weights=[1]):
#only double prec parameters are fit parameters, spec_expt should be len(bgrid)
    xx = []
    for i in range(len(b1_list)):
        b1 = b1_list[i]
        params_in={**params_nonfit,**{x:params_fit[x].value for x in params_fit}}
        params_in['b1'] = b1
        xx.append(cw_spec(bgrid=bgrid[i], params_in=params_in, basis_file='xoxo', prune_on=False)[1]-spec_expt[i])
    return np.hstack(xx) #eps_data

params = Parameters()
scale_init = 1494*1.4/0.9 #366.5 * 0.25 * 10 * 1.05/0.7
shiftg_init = 0
t1edi_init = 5.077 #4.8 + 0.7 -0.3 - 0.3 + 0.15 - 0.07 + 0.1 + 0.03
gib0_init = 1.5
dx_init = 7.5
#shiftx_init = 0

#params.add('b1', value=b1_init, min=0.0005, max=1)
#params.add('shiftg', value=shiftg_init, min=-15, max=15)
params.add('scale', value=scale_init, min=0.001, max=10000)
#params.add('t1edi', value=t1edi_init, min=3, max=8)
#params.add('shiftx', value=shiftx_init, min=0, max=10)
params.add('gib0', value=gib0_init, min=0.01, max=3)
params.add('dx', value=dx_init, min=6, max=9)
params.add('dy', expr='dx')
params.add('dz', expr='dx')
B1max=0.97

dB_list = ['30']
num_spec = len(dB_list)
bgrid = []
spec_expt = []
b1_list = []
weights = []
for db in dB_list:
    #f = 'June28_Saturation_deox/PC5_T39_dB'+db+'.dat'
    f = 'BorisY2021/Apr26/Tempo_DBPh_208K_'+db+'dB.dat'
    aa = np.loadtxt(f,delimiter=',')
    aa = aa[0:-1:8,:]
    weights.append(np.max(aa[:,1])-np.min(aa[:,1]))
    bgrid.append(aa[:,0])
    spec_expt.append(aa[:,1])
    b1_list.append(B1max*10**(-0.05*int(db)))

#weights = [1,1,1,1,1,1]
fix_args={'shiftg':shiftg_init, 't1edi':t1edi_init, 't2edi':0.0, 'nort':0, 'gxx':2.010103, 'gyy':2.0064, 'gzz':2.00233, 'axx':6.11, 'ayy':6.11, 'azz':34.7, 'b0':3320}
out = minimize(sat_residual, params, args=(fix_args, bgrid, spec_expt, b1_list, weights))
report_fit(out)
print('Time taken: ',time.time()-stt)
'''

#boris no potential simulate
'''
d_val = 7.4829
stt=time.time()
scale_init = 2184.507 # 366.5 * 0.25 * 10
shiftg_init = 0
t1edi_init = 5.077#4.80+0.7
B1max=0.97 
dB_list= ['30']
for dbval in dB_list:
    aa=np.loadtxt('BorisY2021/Apr26/Tempo_DBPh_208K_'+dbval+'dB.dat',delimiter=',')
    aa=aa[0:-1:8,:]
    u,v=cw_spec(bgrid=aa[:,0], 
                params_in={'b1':B1max*10**(-0.05*int(dbval)), 
                           'gib0':1.5, 'scale':scale_init, 
                           'shiftg':shiftg_init, 
                           't1edi':t1edi_init, 
                           't2edi':0.0,
                           'dx':d_val,'dy':d_val,'dz':d_val,
                           'nort':0, 'gxx':2.010103, 'gyy':2.0064, 'gzz':2.00233, 'axx':6.11, 'ayy':6.11, 'azz':34.7,
                           'b0':3320
                           },
                basis_file='xoxo', 
                prune_on=0)
    #np.savetxt(file.split('.')[0]+'_guessnort20t1edi5pt077scale1494.txt',np.vstack((u,v)).transpose(), delimiter=',')
    plt.plot(u,v)
    plt.plot(aa[:,0],aa[:,1])
    plt.ylabel('CW spectrum (arb. units)')
    plt.xlabel('Magnetic field (G)')
    plt.legend(['Fit','Actual spectrum'])
    plt.title('$B_1$ = '+np.format_float_positional(B1max*10**(-0.05*int(dbval)), precision=5)+' G, attentuation '+str(dbval)+' dB')
    plt.savefig('TempoDBPh208K_'+file.split('.')[0][10:]+'dB_play.png',fmt='png')
    plt.clf()
print('Time taken: ',time.time()-stt)
'''



# PAPER 3 CALCULATIONS AND PLOTS BEGIN HERE !!!!!!!!!!!
#------------------------------------------
#------------------------------------------

# vs b1
'''
fig, axs = plt.subplots(2, 2)
fig.tight_layout(pad=2.0)
#fig.set_size_inches(20,6)
cnt = 0
for b1 in [0.01, 0.1, 1, 5]:
    u,v=cw_spec(params_in={'gxx':2.0087,
                           'gyy':2.0057,
                           'gzz':2.0021,
                           'axx':6.0,
                           'ayy':6.0,
                           'azz':36.0,
                           'dx':8.0,
                           'dy':8.0,
                           'dz':8.0,
                           'b0': 3360,
                           'b1':b1, 
                           'nort':0,
                           'scale':1,
                           'shiftg':0,
                           't1edi':np.log10(2e5), 
                           't2edi':0.0,
                           'c20':0.0,
                           'oss':0.0
                           },
                basis_file='xoxo', 
                prune_on=0)
    axs[cnt//2, cnt%2].plot(u, v)
    axs[cnt//2, cnt%2].set_title('$B_1$ = '+np.format_float_positional(b1, precision=4)+' G')
    cnt += 1
axs[0,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_xlabel('Magnetic field (G)')
axs[1,1].set_xlabel('Magnetic field (G)')
'''

# pp vs b1
'''
b1_list = []
pp_list = []
for b1 in np.linspace(0.001,2,25):
    u,v=cw_spec(params_in={'gxx':2.0087,
                           'gyy':2.0057,
                           'gzz':2.0021,
                           'axx':6.0,
                           'ayy':6.0,
                           'azz':36.0,
                           'dx':8,
                           'dy':8,
                           'dz':8,
                           'b0': 3360,
                           'b1':b1, 
                           'nort':0,
                           'scale':1,
                           'shiftg':0,
                           't1edi':np.log10(2e5), 
                           't2edi':0.0,
                           'c20':0.0,
                           'oss':0.0
                           },
                basis_file='xoxo', 
                prune_on=0)
    b1_list.append(b1)
    pp_list.append(np.max(v[abs(u-3360)<=10])-np.min(v[abs(u-3360)<=10]))
'''
    
'''
p = Parameters()
p.add('I0', min=0)
p.add('P', min=0)
p.add('eps', min=0, max=10, value=1.5)
def pp_res(params):
    a = np.array(b1_list)
    return np.array(pp_list) - params['I0'] * a/ (1+params['P']*a*a)**params['eps']
out = minimize(pp_res, p)
report_fit(out)
'''


## vs isodiff tensor
'''
fig, axs = plt.subplots(2, 2)
fig.tight_layout(pad=1.0)
#fig.set_size_inches(20,6)
cnt = 0
b1_set = 0.5
for rbar in [7,7.5,8,9]:
    u,v=cw_spec(params_in={'gxx':2.0087,
                           'gyy':2.0057,
                           'gzz':2.0021,
                           'axx':6.0,
                           'ayy':6.0,
                           'azz':36.0,
                           'dx':rbar,
                           'dy':rbar,
                           'dz':rbar,
                           'b0': 3360,
                           'b1':b1_set, 
                           'nort':0,
                           'scale':1,
                           'shiftg':0,
                           't1edi':np.log10(2e5), 
                           't2edi':0.0,
                           'c20':0.0,
                           'oss':0.0
                           },
                basis_file='xoxo', 
                prune_on=0,
                bgrid=np.linspace(-60,60,500)+3360)
    axs[cnt//2, cnt%2].plot(u, v)
    axs[cnt//2, cnt%2].set_title('$R_{xx}=R_{yy}=R_{zz} = 10^{'+str(rbar)+'}\ s^{-1}$')
    cnt += 1
axs[0,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_xlabel('Magnetic field (G)')
axs[1,1].set_xlabel('Magnetic field (G)')
''' 

## vs oss
'''
fig, axs = plt.subplots(2, 2)
fig.tight_layout(pad=1.0)
#fig.set_size_inches(20,6)
cnt = 0
b1_set = 0.5
for oss in [0,7,7.5,8]:
    u,v=cw_spec(params_in={'gxx':2.0087,
                           'gyy':2.0057,
                           'gzz':2.0021,
                           'axx':6.0,
                           'ayy':6.0,
                           'azz':36.0,
                           'dx':8.5,
                           'dy':8.5,
                           'dz':8.5,
                           'b0': 3360,
                           'b1':b1_set, 
                           'nort':0,
                           'scale':1,
                           'shiftg':0,
                           't1edi':np.log10(2e5), 
                           't2edi':0.0,
                           'c20':0.0,
                           'oss':oss
                           },
                basis_file='xoxo', 
                prune_on=0,
                bgrid=np.linspace(-60,60,128)+3360)
    axs[cnt//2, cnt%2].plot(u, v)
    if abs(oss) < 1e-5:
        axs[cnt//2, cnt%2].set_title('$\omega_{ss}= 0$')
    else:
        axs[cnt//2, cnt%2].set_title('$\omega_{ss}= 10^{'+str(oss)+'}\ s^{-1}$')
    cnt += 1
axs[0,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_xlabel('Magnetic field (G)')
axs[1,1].set_xlabel('Magnetic field (G)')
''' 

# dummy fit
'''
bgrid = np.linspace(-60,60,256)+3360
def dummy_sat_residual(params_fit=dict(), params_nonfit=dict(), bgrid=bgrid, spec_expt=np.zeros(128)):
    params_in={**params_nonfit,**{x:params_fit[x].value for x in params_fit}}
    spec_simulated = cw_spec(bgrid=bgrid, params_in=params_in, basis_file='xoxo', prune_on=False)[1]
    return spec_expt-spec_simulated #eps_data
dx_init = 9.5
dz_init = 9.5
params = Parameters()
params.add('scale', value=0.1, min=0.01, max=1)
#params.add('b0', value=2000, min=2000, max=4000)
params.add('dx', value=dx_init, min=5.1, max=10 ,vary=True)
params.add('dy', expr='dx')
params.add('dz', value=dz_init, min=5.1, max=10)
fix_dict = {'nort':0, 'b1':0.5, 'c20':0.0, 'c20':0.0}

# actual model spectrum
xx = cw_spec(bgrid=bgrid, params_in={**fix_dict,**{'scale':1.0, 'dx':7.2, 'dy':7.2, 'dz':8}}, basis_file='xoxo', prune_on=0)[1] #"expt" spectrum

#do lmfit
out = minimize(dummy_sat_residual, params, args=(fix_dict, bgrid, xx), max_nfev=1000)
report_fit(out)


fig, axs = plt.subplots(1, 2)
fig.tight_layout(pad=1.0)
#before fit
axs[0].plot(bgrid, xx)
u,v=cw_spec(bgrid=bgrid, params_in={**fix_dict,**{'scale':0.25, 'dx':dx_init, 'dy':dx_init, 'dz':dz_init}}, basis_file='xoxo', prune_on=0)
ymx = max(max(v),max(xx))
ymn = min(min(v),min(xx))
axs[0].plot(u,v)
axs[0].set_ylabel('CW spectrum (arb. units)')
axs[0].set_xlabel('Magnetic field (G)')
axs[0].set_title('Before least squares fitting')
axs[0].legend(['Model spectrum', 'Initial guess prior to fit'])
#after fit
axs[1].plot(bgrid,xx)
u,v=cw_spec(bgrid=bgrid, params_in={**fix_dict,**{'scale':out.params['scale'].value, 'dx':out.params['dx'].value,
                                                  'dy':out.params['dy'].value, 'dz':out.params['dz'].value}}, basis_file='xoxo', prune_on=0)
ymx = max(max(v),ymx)
ymn = min(min(v),ymn)
axs[1].plot(u,v)
axs[1].set_xlabel('Magnetic field (G)')
axs[1].set_title('After least squares fitting')
axs[1].legend(['Model spectrum', 'Best fit'])

axs[0].set_ylim([ymn,ymx])
axs[1].set_ylim([ymn,ymx])
'''
# just plot the 5PC fits, no fitting  
#'''
fig, axs = plt.subplots(3, 2)
fig.tight_layout(pad=1.0)
stt=time.time()
B1max=0.9 
dB_list= ['30','13','6','2','0']
cnt = 0
for dbval in dB_list:
    aa=np.loadtxt('June28_Saturation_deox/PC5_T39_dB'+dbval+'.dat',delimiter=',')
    aa=aa[0:-1:8,:]
    u,v=cw_spec(bgrid=np.linspace(-150,150,256)+3360,#aa[:,0], 
#scale:  1865.94534 +/- 47.9072161 (2.57%) (init = 1888)
#    t1edi:  5.04481806 +/- 0.01023307 (0.20%) (init = 5.046)
#    gib0:   1.29573435 +/- 0.29837452 (23.03%) (init = 1.94)
#    gib2:   0.65236874 +/- 0.29276701 (44.88%) (init = 0.01)                              
                params_in={'b1':B1max*10**(-0.05*int(dbval)), 
                           'gib0':1.942,#1.5,
                           #'gib2':0.65236874,
                           'shiftg':-3.9, 
                           't1edi':5.046,#5.077, 
                           't2edi':0.0,
                           'nort':20,
                           'c20':2.0,
                           'scale':1887.95325
                           },
                basis_file='xoxo', 
                prune_on=0)
    #np.savetxt(file.split('.')[0]+'_guessnort20t1edi5pt077scale1494.txt',np.vstack((u,v)).transpose(), delimiter=',')
    axs[cnt//2,cnt%2].plot(u,v)
    axs[cnt//2,cnt%2].plot(aa[:,0],aa[:,1])
    axs[cnt//2,cnt%2].legend(['Fit','Actual spectrum'])
    axs[cnt//2,cnt%2].set_title('$B_1$ = '+np.format_float_positional(B1max*10**(-0.05*int(dbval)), precision=5)+' G, attentuation '+str(dbval)+' dB')
    cnt += 1
    #plt.savefig('June28_Saturation_deox/PC5_T39_dB'+dbval+'.pdf',fmt='pdf')
    #plt.clf()
axs[0,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_ylabel('CW spectrum (arb. units)')
axs[1,0].set_xlabel('Magnetic field (G)')
axs[1,1].set_xlabel('Magnetic field (G)')    
    
axs[2,0].set_xlabel('Magnetic field (G)')    
axs[2,0].set_ylabel('CW spectrum (arb. units)')
print('Time taken: ',time.time()-stt)
#plt.show()
#'''

#fit Boris 5PC data
'''
stt=time.time()
def sat_residual(params_fit=dict(), params_nonfit=dict(), bgrid=np.reshape(np.linspace(-60,60,256)+3360,(1,-1)), spec_expt=np.zeros((1,128)), b1_list = [0.1], weights=[1]):
#only double prec parameters are fit parameters, spec_expt should be len(bgrid)
    xx = []
    for i in range(len(b1_list)):
        b1 = b1_list[i]
        params_in={**params_nonfit,**{x:params_fit[x].value for x in params_fit}}
        params_in['b1'] = b1
        xx.append(weights[i]*(cw_spec(bgrid=bgrid[i], params_in=params_in, basis_file='xoxo', prune_on=False)[1]-spec_expt[i]))
    return np.hstack(xx) #eps_data

params = Parameters()
scale_init = 1888 #1494 #366.5 * 0.25 * 10 * 1.05/0.7
shiftg_init = -3.9
t1edi_init = 5.046 #5.077 #4.8 + 0.7 -0.3 - 0.3 + 0.15 - 0.07 + 0.1 + 0.03
gib0_init = 1.94#1.5
gib2_init = 0.01
#shiftx_init = 0

#params.add('b1', value=b1_init, min=0.0005, max=1)
#params.add('shiftg', value=shiftg_init, min=-15, max=15)
params.add('scale', value=scale_init, min=0.001, max=10000)
params.add('t1edi', value=t1edi_init, min=3, max=8)
#params.add('shiftx', value=shiftx_init, min=0, max=10)
params.add('gib0', value=gib0_init, min=0.01, max=3)
params.add('gib2', value=gib2_init, min=0, max=3)
B1max=0.9

dB_list = ['30','13','6','2','0'] #['0','2','4','10','20','40']
num_spec = len(dB_list)
bgrid = []
spec_expt = []
b1_list = []
weights = []

for db in dB_list:
    f = 'June28_Saturation_deox/PC5_T39_dB'+db+'.dat'
    aa = np.loadtxt(f,delimiter=',')
    aa = aa[0:-1:8,:]
    bgrid.append(aa[:,0])
    spec_expt.append(aa[:,1])
    b1_list.append(B1max*10**(-0.05*int(db)))
    weights.append(1/(max(aa[:,1])-min(aa[:,1])))

out = minimize(sat_residual, params, args=({'shiftg':shiftg_init, 'nort':20}, bgrid, spec_expt, b1_list, weights))
report_fit(out)
print('Time taken: ',time.time()-stt)
'''

'''
FIT AGAIN ON APR 27, TOOK ~11k seconds, 10 orientations
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 21
    # data points      = 640
    # variables        = 3
    chi-square         = 0.90174186
    reduced chi-square = 0.00141561
    Akaike info crit   = -4195.53290
    Bayesian info crit = -4182.14850
[[Variables]]
    scale:  1887.95325 +/- 48.6151116 (2.58%) (init = 1494)
    t1edi:  5.04603352 +/- 0.01036494 (0.21%) (init = 5.077)
    gib0:   1.94223282 +/- 0.08735447 (4.50%) (init = 1.5)
[[Correlations]] (unreported correlations are < 0.100)
    C(scale, t1edi) = -0.842
    C(scale, gib0)  =  0.440
'''

''' 
FIT AGAIN ON APR 25, TOOK 13332 seconds, 8 orientations
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 25
    # data points      = 768
    # variables        = 3
    chi-square         = 107.610454
    reduced chi-square = 0.14066726
    Akaike info crit   = -1503.32884
    Bayesian info crit = -1489.39747
[[Variables]]
    scale:  1790.24558 +/- 36.4245140 (2.03%) (init = 1494)
    t1edi:  5.08824387 +/- 0.01085657 (0.21%) (init = 5.077)
    gib0:   1.53034846 +/- 0.06546246 (4.28%) (init = 0.5)
[[Correlations]] (unreported correlations are < 0.100)
    C(scale, t1edi) = -0.802
    C(scale, gib0)  =  0.404
Time taken:  13332.69223189354
'''

'''
FIT AGAIN ON APR 28, TOOK ~55k seconds (nort=20 instead of 10), 4 parameters (gib2 added)
[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 53
    # data points      = 640
    # variables        = 4
    chi-square         = 0.88989694
    reduced chi-square = 0.00139921
    Akaike info crit   = -4201.99539
    Bayesian info crit = -4184.14952
[[Variables]]
    scale:  1865.94534 +/- 47.9072161 (2.57%) (init = 1888)
    t1edi:  5.04481806 +/- 0.01023307 (0.20%) (init = 5.046)
    gib0:   1.29573435 +/- 0.29837452 (23.03%) (init = 1.94)
    gib2:   0.65236874 +/- 0.29276701 (44.88%) (init = 0.01)
[[Correlations]] (unreported correlations are < 0.100)
    C(gib0, gib2)   = -0.956
    C(scale, t1edi) = -0.829
    C(scale, gib0)  =  0.281
    C(scale, gib2)  = -0.163
Time taken:  55056.00435447693
'''