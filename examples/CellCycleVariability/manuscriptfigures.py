import numpy as np
import seaborn as sns
from scipy.stats import poisson,nbinom,binom,skew,kurtosis
import matplotlib.pyplot as plt
import runcyclo as gill
import erlanggenerating as erl
from scipy.special import factorial 
from scipy.optimize import minimize
from matplotlib.ticker import ScalarFormatter, NullFormatter
from tqdm import tqdm
import bursttrans as bt
import sympy as sp
import matplotlib.gridspec as gridspec # GRIDSPEC !
from matplotlib.colors import LogNorm
from scipy.interpolate import UnivariateSpline
import pandas as pd

sbcolorcyclelight = sns.color_palette('deep')
sbcolorcycledark = sns.color_palette('dark')
sbcolorcycleblue = sns.color_palette('GnBu_d')

plt.ioff()

#####
# Dictionary with the parameters of the model
pardict = {
"Omega": 1,
"k0" : 10,  "k1" : 1, "k2" : 0, "k3" : 0, "T" : 1,
"m0" : 0, "p0" : 0, "ffwrite" : 1, "timelapse": 0.01,
"totaltime": 1000000, "dt" : 0.01, "runtype" : 0, "SEED": -1,
"runtimes": 1, "stocycle":1, "phasenumber":2, "presynthesis":1,
"cellphaserates" :[1.0] 
}


############################################################################
###### The theoretical expressions are calculated with the implementation 
# in this script or with the fastest code of burstrans.py
# the stochast simulations are managed with gillespie_cyclo.py

#################################################################################
#
# Theoretical expressions for the mean number of mRNA in the Erlang distribution
# for trajectory and population measurements
# star denotes deterministic counterparts
#
def  meantrajErlang(w,nhat,eta,Delta): # theoretical expression

    factor = 1/(1+eta*Delta)
    factor1 = factor**((1-w)/Delta)
    factor2 = factor**(1/Delta)
    avern = 1-factor1/(2-factor2) 
    avern *= -1*nhat/eta
    avern += (1-w)*2*nhat
    avern += w*nhat
    return avern

def meanstartrajErlang(w,nhat,eta,Delta):
    factor = np.exp(-eta*(1-w))
    factor /= 2-np.exp(-eta)
    avern = w*nhat + (1-w)*2*nhat-nhat/eta*(1-factor)
    return avern

def meanpopErlang(w,nhat,eta,Delta):
    factor = 2**(1-w)*Delta*eta*nhat
    factor /= (-1+(2**Delta)+Delta*eta)
    return factor

def meanstarpopErlang(w,nhat,eta,Delta):
    factor = 2**(1-w)*nhat*eta
    factor /= (eta+np.log(2))
    return factor
#####################################################################################


#####################################################################################
# 
# General expressions for the factorial moments and statistics of the hypoexponential distribution
# for faster results most of the routines use the ones defined in burstrans.py


def factmom_Erlang(j,p,nhat,kd,W,N): # (n_j)_p
    # index j of delta for kd = k/d, nhat = r_0/d and factorial moment p for Erlang

    if p==0:
        return 1/N
    elif j==1:
        numerator = 2*nhat*p*factmom_Erlang(1,p-1,nhat,kd,W,N) + kd/2**(p-1)*thetaj(N,p,nhat,kd,W,N)
        denominator = 2*(p+kd)-kd/2**(p-1)*deltaj(N,p,kd)
        return numerator/denominator
    else:
        return deltaj(j,p,kd)*factmom_Erlang(1,p,nhat,kd,W,N)+thetaj(j,p,nhat,kd,W,N)

def deltaj(j,p,kd):
    # index j of delta for kd = k/d and factorial moment p for Erlang
    return (kd/(kd+p))**(j-1)   


def thetaj(j,p,nhat,kd,W,N):
    # index j of theta for kd = k/d and nhat = r(0)/d for factorial moment p for Erlang
    sum1 = 0
    for m in range (1,j): # until j-1
        if m<W:
            rfactor = 1
        else:
            rfactor = 2
        sum1 += rfactor*factmom_Erlang(m+1,p-1,nhat,kd,W,N)/deltaj(m+1,p,kd)

    return deltaj(j,p,kd)*nhat*p/(kd+p) * sum1

def meanErlang(nhat,eta,W,N,mode='traj'):

    if mode == 'traj':
        pi = 1/N*np.ones(N) # probability of age i 
        k = N*np.ones(N)
        kd = k/eta # k/d (assumes T=1)
    elif mode == 'pop':
        pi = (2**(1/N)-1)/2**(np.arange(1,N+1)/N-1)
        k = N*np.ones(N) * 2**(1/N)
        kd = k/eta
    else: 
        print('Unknown mode of observation')

    sum1 = 0
    for m in np.arange(1,N+1):
        sum1+= pi[m-1]*factmom_Erlang(m,1,nhat,kd[m-1],W,N)/factmom_Erlang(m,0,nhat,kd[m-1],W,N)

    return sum1

def varstartrajErlang(nhat,eta,w,Delta):
    # corresponding with "varav" in Ramon's Mathematica's notebook

    varstar = 3*nhat*w + w - 4*nhat -2
    varstar = varstar*2*eta + 7*nhat + 2
    varstar = varstar*2*np.exp(eta) - (3*nhat*w + w - 4*nhat -2)*2*eta
    varstar = -2-5*nhat + varstar
    varstar = varstar*(-1+2*np.exp(eta))
    varstar = -3*np.exp(2*eta*w)*nhat + 2*np.exp(eta*w)*(-1+2*np.exp(eta))*(1+6*nhat)-varstar
    varstar *= 1/((1-2*np.exp(eta))**2*eta)
    varstar += -2*nhat*(-2+1/eta+np.exp(eta*w)/(eta-2*np.exp(eta)*eta)+w)**2
    varstar *= 0.5*nhat

    return varstar

####################################################################################################


def getFig23(panel='a',figure=2,loaddata=False):
   # Code to obtain figures 2 and 3 of the manuscript 
   # figure = 2 contains trajectory measurements
   # figure = 3 contains population measurements

   # panel names correspond with the order of figure 1
    
    if panel=='a':
        conditionlist = [{'N':2,'W':1,'nhat':1,'label':'$\\hat{n}=1$','mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':500,'label':'$\\hat{n}=500$', 'mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':500,'label':'$\\hat{n}=500$', 'mRNAgeo':10}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)


    if panel=='b':
        conditionlist = [{'N':12,'W':3,'nhat':50,'label':'$\\Delta=1/12\\quad w=1/4$'},{'N':12,'W':6,'nhat':50,'label':'$\\Delta=1/12\\quad w=1/2$'},{'N':12,'W':9,'nhat':50,'label':'$\\Delta=3/4\\quad w=3/4$'}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)

    if panel=='d':
        conditionlist = [{'N':2,'W':1,'nhat':1,'label':'$\\hat{n}=1$','mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':500,'label':'$\\hat{n}=500$', 'mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':500,'label':'$\\hat{n}=500$', 'mRNAgeo':10}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)

    if panel=='e':
        conditionlist = [{'N':12,'W':3,'nhat':50,'label':'$\\Delta=1/12\\quad w=1/4$','mRNAgeo':-1},
                        {'N':12,'W':6,'nhat':50,'label':'$\\Delta=1/12\\quad w=1/2$','mRNAgeo':-1},
                        {'N':12,'W':9,'nhat':50,'label':'$\\Delta=1/12\\quad w=3/4$','mRNAgeo':-1}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)

    if panel=='f': # new panel for figure 2 showing differences in error depending on burst size
        conditionlist = [{'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':1},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':10},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':100},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$','mRNAgeo':-1}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)


    simulationtime = 500
    simulationburnout =  100
    if (figure in [2] and panel in ['a','b','f']):
        repeats = 50
    elif (figure in [2] and panel in ['d','e']):
        repeats = 200
    elif (figure in [3] and panel in ['a','c']):
        repeats = 5000
    elif (figure in [3] and panel in ['b','e']):
        repeats = 25000

    for icondition,condition in enumerate(conditionlist):

        w = condition['W']/condition['N']
        Delta = 1/condition['N']

        etanum_corrected = etanum*(1+0.15*(icondition-1))

        if (condition['mRNAgeo']== -1 and figure==2 and panel in ['d','e','f']): # variance in constitutive expression in a trajectory
            mstartheo = meanstartrajErlang(w,1,etatheo,Delta)
            mtheo = meantrajErlang(w,1,etatheo,Delta)
            R = (mtheo-mstartheo)/mtheo
            Ndet = 400
            k = condition['N']
            varstartheo = [bt.vartrajErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                            Ndet,etatheo_i,int(condition['W']/condition['N']*Ndet),Ndet,
                                            mode = 'constitutive') for etatheo_i in etatheo]
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            vartheo = [bt.vartrajErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                           k,etatheo_i,condition['W'],condition['N'],
                                           mode = 'constitutive') for etatheo_i in etatheo]
            vartheo = np.array(vartheo)
            varstartheo = np.array(varstartheo)
            varR = (vartheo-varstartheo)/vartheo
            print('varstarteo',varstartheo)
            print('vartheo',vartheo)
            print('varR',varR)

        elif (condition['mRNAgeo'] > 0 and figure==2 and panel in ['d','e','f']): # variance in bursty expression in a trajectory
            mstartheo = meanstartrajErlang(w,1,etatheo,Delta)
            mtheo = meantrajErlang(w,1,etatheo,Delta)
            R = (mtheo-mstartheo)/mtheo
            k = condition['N']
            Ndet = 400
            varstartheo = [bt.vartrajErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                            Ndet,etatheo_i,int(condition['W']/condition['N']*Ndet),Ndet) for etatheo_i in etatheo]
            varstartheo = np.array(varstartheo)
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            vartheo = [bt.vartrajErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                           k,etatheo_i,condition['W'],condition['N']) for etatheo_i in etatheo]
            vartheo = np.array(vartheo)
            varR = (vartheo-varstartheo)/vartheo
            print('varstarteo',varstartheo)
            print('vartheo',vartheo)
            print('varR',varR)

        elif (condition['mRNAgeo']== -1 and figure==3  and panel in ['a','b']): # mean constitutive expression in a population
            
            Ndet = 200
            k = condition['N']
            meanstartheo = [bt.meanpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                            Ndet,etatheo_i,int(condition['W']/condition['N']*Ndet),Ndet,
                                            mode = 'constitutive') for etatheo_i in etatheo]
            meanstartheo = np.array(meanstartheo)
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            meantheo = [bt.meanpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                           k,etatheo_i,condition['W'],condition['N'],
                                           mode = 'constitutive') for etatheo_i in etatheo]
            meantheo = np.array(meantheo)
            meanR = (meantheo-meanstartheo)/meantheo
            print('varstarteo',meanstartheo)
            print('vartheo',meantheo)
            print('varR',meanR)

        elif (condition['mRNAgeo']== -1 and figure==3 and panel in ['d','e','f']): # variance constitutive expression in a population

            Ndet = 1000
            k = condition['N']
            varstartheo = [bt.varpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                            Ndet,etatheo_i,
                                            int(condition['W']/condition['N']*Ndet),Ndet,
                                            mode = 'constitutive') for etatheo_i in etatheo]
            varstartheo = np.array(varstartheo)
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            vartheo = [bt.varpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                           k,etatheo_i,condition['W'],condition['N'],
                                           mode = 'constitutive') for etatheo_i in etatheo]
            vartheo = np.array(vartheo)
            varR = (vartheo-varstartheo)/vartheo
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            print('varstarteo',varstartheo)
            print('varstar',vartheo)
            print('varR',varR)

        elif (condition['mRNAgeo']>0 and figure==3 and panel in ['d','e','f']): # variance bursty expression in a population
            
            Ndet = 200
            k = condition['N']
            varstartheo = [bt.varpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                            Ndet,etatheo_i,int(condition['W']/condition['N']*Ndet),Ndet) for etatheo_i in etatheo]
            varstartheo = np.array(varstartheo)
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            vartheo = [bt.varpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                           k,etatheo_i,condition['W'],condition['N']) for etatheo_i in etatheo]
            vartheo = np.array(vartheo)
            varR = (vartheo-varstartheo)/vartheo
            print('varstarteo',varstartheo)
            print('vartheo',vartheo)
            print('varR',varR)

        elif (condition['mRNAgeo']>0 and figure==3 and panel in ['a','b']): # mean burst expression population
            
            Ndet = 200
            k = condition['N']
            meanstartheo = [bt.meanpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                            Ndet,etatheo_i,int(condition['W']/condition['N']*Ndet),Ndet) for etatheo_i in etatheo]
            varstartheo = np.array(varstartheo)
            print(condition['nhat'],etatheo,w,condition['W'],Delta)
            meantheo = [bt.meanpopErlang(condition['mRNAgeo'],condition['nhat']*etatheo_i,
                                           k,etatheo_i,condition['W'],condition['N']) for etatheo_i in etatheo]
            meantheo = np.array(meantheo)
            meanR = (meantheo-meanstartheo)/meantheo
            print('varstarteo',meanstartheo)
            print('vartheo',meantheo)
            print('varR',meanR)

        # Plotting theoretical prediction for R
        if (panel in ['a','b']):
                plt.plot(etatheo,meanR,'-',color = sbcolorcyclelight[icondition],label='theory')
        else:
            plt.plot(etatheo,varR,'-',color = sbcolorcyclelight[icondition], label ='theory')

        savelist = []

        # Simulation for constitutive expression
        if loaddata: # if loaddata is true and prior simulations exist, they are recycled
            D = np.loadtxt('fig{}_panel{}_cond{}.dat'.format(figure,panel,icondition))
        else:
            for eta in etanum_corrected:
                    # simulation time
                    pardict['phasenumber'] = condition['N']
                    pardict['presynthesis'] = condition['W']
                    pardict['stocycle'] = 1
                    pardict['mRNAgeo'] = condition['mRNAgeo']
                    pardict['k1'] = eta # because cell cucle lenght is T =1
                    pardict['k0'] = condition['nhat']*eta/np.abs(condition['mRNAgeo']) # abs keeps constitutive =1
                    pardict['m0'] = condition['nhat']
                    tau = 1 # total cell cycle duration
                    pardict["T"] = tau
                    pardict["cellphaserates"] = np.ones(condition['N'])/(tau/condition['N'])
                    print("Calculating condition eta:{}, condition:{}".format(eta,icondition))

                    marray = np.ones(0)         
                    mvararray = np.ones(0)          

                    for repeat in tqdm(range(repeats)):

                    # Simulations for the stochastic cell cycle case
                        if figure in [2]:
                            simulationtime = 600*max(1/pardict['k1'],1/pardict['k0'],tau)
                            simulationburnout = 0.4*simulationtime
                            pardict['totaltime'] = simulationtime
                            pardict['ffwrite'] = 0
                            timelapse = simulationtime/100000
                            pardict['timelapse'] = timelapse
                            B = gill.RunTime_C(pardict)
                            sliceburnout = B[:,0]>simulationburnout
                            m,p = averagefromtraj_timelapse(B[sliceburnout])
                            mvar,pvar = varfromtraj(B[sliceburnout])
                        elif figure in [3]:
                            pardict['ffwrite'] = 4 # only write last point (that will onw per cell)
                            simulationtime = 10
                            pardict['totaltime'] = simulationtime
                            pardict['m0'] = meanErlang(condition['nhat'],eta,condition['W'],condition['N'],'pop')
                            B = gill.RunTime_C_Pop(pardict)
                            m,p = np.mean(B[:,1]),np.mean(B[:,2])
                            mvar,pvar = np.var(B[:,1]),np.var(B[:,2])

                        marray = np.append(marray,m)
                        mvararray = np.append(mvararray,mvar)


                    # Simulations for the deterministic cell cycle case
                    if ((w>0) and (w<1)):

                        pardict['stocycle'] = 0
                        pardict['phasenumber'] = 2
                        pardict['presynthesis'] = 1
                        pardict["cellphaserates"] = [1/w,1/(1-w)]

                    elif w==0:

                        pardict['stocycle'] = 0
                        pardict['phasenumber'] = 1
                        pardict['presynthesis'] = 0
                        pardict["cellphaserates"] = [1]

                    else: # w=1

                        pardict['stocycle'] = 0
                        pardict['phasenumber'] = 1
                        pardict['presynthesis'] = 1
                        pardict["cellphaserates"] = [1]

                    mstararray = np.ones(0)
                    mstarvararray = np.ones(0)
                    for repeat in tqdm(range(repeats)):
                        if figure in [2]:
                            simulationtime = 600*max(1/pardict['k1'],1/pardict['k0'],tau)
                            simulationburnout = 0.4*simulationtime
                            pardict['totaltime'] = simulationtime
                            pardict['ffwrite'] = 0
                            timelapse = simulationtime/100000
                            pardict['timelapse'] = timelapse
                            B = gill.RunTime_C(pardict)
                            sliceburnout = B[:,0]>simulationburnout
                            m,p = averagefromtraj_timelapse(B[sliceburnout])
                            mvar,pvar = varfromtraj(B[sliceburnout])
                        elif figure in [3]:
                            pardict['ffwrite'] = 4 # only write last point (that will onw per cell)
                            urnd = np.random.random()
                            #simulationtime = 10+np.random.random()*pardict["T"] # to ensure that there is no synchornization between simulations for the deterministic case
                            simulationtime = 10 + np.log(2/(2-urnd))/np.log(2) # random time initating at a random cell of a population
                            pardict['totaltime'] = simulationtime
                            pardict['m0'] = meanErlang(condition['nhat'],eta,condition['W'],condition['N'],'pop')
                            B = gill.RunTime_C_Pop(pardict)
                            m,p = np.mean(B[:,1]),np.mean(B[:,2])
                            mvar,pvar = np.var(B[:,1]),np.var(B[:,2])
                        mstararray = np.append(mstararray,m)
                        mstarvararray = np.append(mstarvararray,mvar)

                    meanmRNA = np.mean(marray) # mean of the mean of the population
                    varmRNA = np.mean(mvararray) # mean of the var of population

                    varmeanmRNA = np.var(marray) # var of the mean of the population
                    varvarmRNA = np.var(mvararray) # var of the var of the population

                    meanmRNAstar = np.mean(mstararray) # mean of the mean of the population
                    varmRNAstar = np.mean(mstarvararray) # mean of the var of population

                    varmeanmRNAstar = np.var(mstararray) # var of the mean of the population
                    varvarmRNAstar = np.var(mstarvararray) # var of the var of the population

                    savelist.append([eta,repeats,
                        meanmRNA,varmRNA,varmeanmRNA,varvarmRNA,
                        meanmRNAstar,varmRNAstar,varmeanmRNAstar,varvarmRNAstar])

            np.savetxt('fig{}_panel{}_cond{}.dat'.format(figure,panel,icondition),np.array(savelist))
            D = np.array(savelist)

        etanum_corrected = D[:,0]
        replicates = D[:,1]

        meanmRNA = D[:,2]
        varmRNA = D[:,3]
        varmeanmRNA = D[:,4] 
        varvarmRNA = D[:,5] 

        meanmRNAstar = D[:,6] 
        varmRNAstar = D[:,7]
        varmeanmRNAstar = D[:,8]
        varvarmRNAstar = D[:,9]

        if panel in ['a','b']:

            rep = replicates

            meanR = (meanmRNA- meanmRNAstar)/(meanmRNA)
            err_relmeanmRNA = np.sqrt(varmeanmRNA/rep)/meanmRNA
            err_relmeanmRNAstar = np.sqrt(varmeanmRNAstar/rep)/meanmRNAstar

            # errR = (1+meanR)/np.sqrt(2*(replicates-1)*10)
            errR = (meanR-1)*np.sqrt(err_relmeanmRNA**2 + err_relmeanmRNAstar**2)

            plt.plot(etanum_corrected,meanR,'o',color = (0,0,0,0),
                markeredgewidth = 1.5,markersize = 4, markeredgecolor = sbcolorcyclelight[icondition],label=condition['label'])
            plt.errorbar(etanum_corrected,meanR,yerr =1.5*errR,ls='none',color=sbcolorcyclelight[icondition])#fmt=None,color = sbcolorcycleblue[icondition])

            # plt.plot(etanum_corrected,meanmRNA,'o',color = (0,0,0,0),
            #   markeredgewidth = 1.5,markersize = 4, markeredgecolor = sbcolorcyclelight[icondition],label=condition['label'])
            # plt.errorbar(etanum_corrected,meanmRNA,yerr =1.5*errR,ls='none',color=sbcolorcyclelight[icondition])#fmt=None,color = sbcolorcycleblue[icondition])


        else:

            Rvar = (varmRNA - varmRNAstar)/varmRNA

            varmRNA_corrected = varmRNA + varmeanmRNA
            varmRNAstar_corrected = varmRNAstar + varmeanmRNAstar
            Rvar_corrected = (varmRNA_corrected - varmRNAstar_corrected)/varmRNA_corrected

            err_relvarmRNA = np.sqrt(varvarmRNA/replicates)/varmRNA_corrected
            err_relvarmRNAstar = np.sqrt(varvarmRNAstar/replicates)/varmRNAstar_corrected

            # errRvar_corrected = (1-Rvar_corrected)*(err_relvarmRNA + err_relvarmRNAstar)
            errRvar_corrected = (1+Rvar_corrected)/np.sqrt(2*(replicates-1))

            plt.plot(etanum_corrected,Rvar_corrected,'o',color = (0,0,0,0),
                markeredgewidth = 1.5,markersize = 4, markeredgecolor = sbcolorcyclelight[icondition],label=condition['label'])
            plt.errorbar(etanum_corrected,Rvar_corrected,yerr =errRvar_corrected,ls='none',color=sbcolorcyclelight[icondition])#fmt=None,color = sbcolorcycleblue[icondition])          


        plt.legend()
    
    plt.xlabel('$\\eta$')
    if panel in ['d','e','f']:
        plt.ylabel('$R_\\sigma$')
        #plt.ylabel('$\\sigma^2$')
    else:
        plt.ylabel('$R$')
    plt.xscale('log')
    ax = plt.gca()
    plt.tight_layout()
    plt.savefig('fig{}{}.pdf'.format(figure,panel))
    plt.show()

def getFig3a(output = 'trajhistogram'):

    if output == 'trajhistogram':
        eta = 1.0
        nhat = 50
        trajrepeats = 100
        pardict['totaltime'] = 8

    elif output == 'stagehistogram':
        eta = 0.00001
        nhat = 1.0  
        trajrepeats = 2000
        pardict['totaltime'] = 12

    N=4
    pardict['k1'] = eta
    pardict['k0'] = nhat*eta
    pardict['m0'] = int(meantrajErlang(0.75,nhat,eta,1/N))
    pardict['ffwrite'] = 1
    # Simulations of population
    pardict['stocycle'] = 1 # pardict['k1'] = 0.000001 # because cell cucle lenght is T =1
    #simulationtime = 600*max(4/pardict['k1'],4/pardict['k0'],tau)
    pardict['presynthesis'] = 2 ## the actual condition is n<presynthesis
    pardict["cellphaserates"] = np.ones(N)/(1/N)
    print("Single traj loop for eta={} and N ={}".format(eta,N))

    mRNA_histo_traj = np.ones(0)
    mRNA_stage_traj = np.ones(0)

    # Simulation for a trajectory of stochastic cell cycle
    for n in tqdm(range(trajrepeats)):
        C = gill.RunTime_C(pardict) 
        mRNA_histo_traj = np.concatenate((mRNA_histo_traj,C[:,1]))
        mRNA_stage_traj = np.concatenate((mRNA_stage_traj,C[:,3]))

    # # Simulation for a trajectory of deterministic cell cycle
    # for n in tqdm(range(trajrepeats)):
    #   pardict['stocycle'] = 0
    #   pardict['phasenumber'] = 4
    #   pardict['presynthesis'] = 2
    #   pardict["cellphaserates"] = [1/2,1/2]
    #   C = gill.RunTime_C(pardict) 
    #   mRNA_histo_traj = np.concatenate((mRNA_histo_traj,C[:,1]))
    #   mRNA_stage_traj = np.concatenate((mRNA_stage_traj,C[:,3]))

    # urnd = np.random.random()
    # simulationtime = 10 + np.log(2/(2-urnd))/np.log(2) # random time initating at a random cell of a population


    if output == 'stagehistogram':
        pardict['ffwrite'] = 4 # Only last snapshot


    C = gill.RunTime_C_Pop(pardict)
    time_pop = C[:,0]
    time_pop = np.concatenate((np.zeros(1),time_pop))
    mRNA_pop = C[:,1]
    mRNA_stage_pop = C[:,3]
    cuts = np.where((time_pop[1:]-time_pop[:-1])<0)[0]

    time_single = C[:cuts[0],0]
    traj_single = C[:cuts[0],1]

    if output == 'trajhistogram':

        fig = plt.figure(1, figsize=(10,8))
        gs = gridspec.GridSpec(1,2, height_ratios=[1], width_ratios=[0.7,0.1])
        gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.03)
        ax_traj = plt.subplot(gs[0,0]) # place it where it should be.
        ax_histo = plt.subplot(gs[0,1]) # place it where it should be.

        print('cuts', cuts)
        for icut,cut in enumerate(cuts[:-1]):
                nextcut = cuts[icut+1]
                print("cutpair", cut, nextcut,time_pop[nextcut-1])
                ax_traj.plot(time_pop[cut+1:nextcut]+(pardict['totaltime']-time_pop[nextcut-1]), # we add the time necessary to set the t=0 of the trajectory with its time in the population division
                            mRNA_pop[cut+1:nextcut],
                            '-',color = sbcolorcyclelight[1], alpha = 0.5)

        ax_traj.plot(time_single,traj_single,'-',color = sbcolorcycledark[0],lw = 1.5)

        ax_traj.set_xlabel('time')
        ax_traj.set_ylabel('mRNA number')

        ax_histo.hist(mRNA_histo_traj,bins = 20, orientation = 'horizontal', density = True, 
            color = sbcolorcycledark[0], alpha = 0.3)
        ax_histo.hist(mRNA_histo_traj,bins = 20, orientation = 'horizontal', histtype = 'step', density = True,
            edgecolor = sbcolorcycledark[0], fc = 'None', lw = 1.5)

        ax_histo.hist(mRNA_pop[cuts-1],bins = 20, orientation = 'horizontal', density = True,
            color = sbcolorcycledark[1], alpha = 0.3)
        ax_histo.hist(mRNA_pop[cuts-1],bins = 20, orientation = 'horizontal', histtype = 'step',density = True,
            edgecolor = sbcolorcycledark[1], fc = 'None', lw = 1.5)

        ax_histo.set_xticklabels([])
        ax_histo.set_yticklabels([])

        ax_traj.set_ylim([0,130])
        ax_histo.set_ylim([0,130])

        plt.savefig('histogram.pdf')

    elif output == 'stagehistogram':


        xb = np.arange(1,5)
        # plotting theo homogeneous
        
        plt.bar(xb - 0.15, np.ones(4)*0.25, width = 0.3, color = sbcolorcycledark[0], alpha = 0.4)

        # plotting theo pop     
        pi = (2**(1/N)-1)/(2**(xb/N-1))
        plt.bar(xb + 0.15, pi, width = 0.3, color = sbcolorcycledark[1], alpha = 0.4)

        num_traj, bins = np.histogram(mRNA_stage_traj+1, bins=[0.5,1.5,2.5,3.5,4.5], density = True)
        num_pop, bins = np.histogram(mRNA_stage_pop+1,   bins=[0.5,1.5,2.5,3.5,4.5], density = True)

        print(num_traj)
        plt.plot(np.arange(1,5)-0.15,num_traj,'o',color = sbcolorcycledark[0], lw = 0)
        plt.plot(np.arange(1,5)+0.15,num_pop,'o',color = sbcolorcycledark[1], lw = 0)
        plt.errorbar(np.arange(1,5)-0.15,num_traj, np.sqrt(num_traj*(1-num_traj)/(trajrepeats*pardict['totaltime']/N)), color = sbcolorcycledark[0], ls = 'none')
        plt.errorbar(np.arange(1,5)+0.15,num_pop, 2.0*np.sqrt(num_pop*(1-num_pop)/len(mRNA_stage_pop)), color = sbcolorcycledark[1], ls = 'none')

        plt.xticks([1,2,3,4],[1,2,3,4])
        plt.yticks([0.1,0.2,0.3],[0.1,0.2,0.3])

        plt.savefig('stages.pdf')


    plt.show()


def getFig1():

    # Figure 1
    fig = plt.figure(1, figsize=(15,8))
    gs = gridspec.GridSpec(3,1)

    gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.00)

    nhat = 50
    # pardict['k0'] = 10.0
    tau = 1 # total cell cycle duration
    pardict["T"] = tau

    Ns = [4,4,4]
    Ws = [3,3,3]
    etas = [0.1,1,10]
    lims = [[0,20],[0,100],[0,150]]

    ax_trajs = [0,0]
    for ieta,eta in enumerate(etas):
        N = Ns[ieta]
        pardict['phasenumber'] = N
        pardict['k1'] = eta
        pardict['k0'] = nhat*eta

        pardict['m0'] = int(meantrajErlang(0.75,nhat,eta,1/N))
        pardict['presynthesis'] = Ws[ieta] ## the actual condition is n<presynthesis

        totaltime = 30

        # Simulations with stochastic cell cycle
        pardict['ffwrite'] = 1
        pardict['stocycle'] = 1 # pardict['k1'] = 0.000001 # because cell cucle lenght is T =1
        #simulationtime = 600*max(4/pardict['k1'],4/pardict['k0'],tau)
        pardict['totaltime'] = totaltime
        pardict["cellphaserates"] = np.ones(N)/(1/N)
        print("Single traj loop for eta={} and N ={}".format(eta,N))
        B = gill.RunTime_C(pardict)
        time_stochastic = B[:,0]
        time_stochastic = np.concatenate((np.zeros(1),time_stochastic))
        cycle_stochastic = B[:,1]
        cycle_stochastic = np.concatenate((np.ones(1)*cycle_stochastic[0],cycle_stochastic))

        # Simulations with fixed cell cycle
        pardict['ffwrite'] = 1
        pardict['stocycle'] = 0 # pardict['k1'] = 0.000001 # because cell cucle lenght is T =1
        #simulationtime = 600*max(4/pardict['k1'],4/pardict['k0'],tau)
        pardict['totaltime'] = totaltime
        pardict["cellphaserates"] = np.ones(N)/(1/N)
        print("Single traj loop for eta={} and N ={}".format(eta,N))
        C = gill.RunTime_C(pardict)
        time_fixed = C[:,0]
        time_fixed = np.concatenate((np.zeros(1),time_fixed))
        cycle_fixed = C[:,1]
        cycle_fixed = np.concatenate((np.ones(1)*cycle_fixed[1],cycle_fixed))

        ax_trajs[ieta] = plt.subplot(gs[ieta,0]) # place it where it should be.
        
        # for i in range(totaltime):
        #     ax_trajs[ieta].plot([i,i],[0,lims[ieta][1]],'k-')
        
        ax_trajs[ieta].step(time_fixed,cycle_fixed,where='post')
        ax_trajs[ieta].step(time_stochastic,cycle_stochastic,where='post')
        ax_trajs[ieta].set_ylim([0,lims[ieta][1]])
        ax_trajs[ieta].set_xlim([0,pardict['totaltime']])

        ax_trajs[ieta].set_xticks([])
        ax_trajs[ieta].set_xticklabels([])

    ax_trajs[2].set_xticks([0.6666,1,2,3,4,5]) 
    ax_trajs[2].set_xticklabels(['W','T','2T','3T','4T','5T'])

    plt.xlabel('time')
    ax_trajs[1].set_ylabel('mRNA number')

    plt.savefig('trajseta_Fig2.pdf')
    plt.show()


def getFig4(mode = 'lineage', burst = -1, output = 'grid', loaddata = False, N = 12, showcomponents = False):

    # to recover the different panels 
    # mode can be  'lineage' or 'population'
    # output can be 'grid' (panel b),  'line' (panel c), or  'distribution' (panel a)



    # conditionlist = [{'eta':0.1,'nhat':50,'N':2},
    #               {'eta':1.0,'nhat':50,'N':2},
    #               {'eta':10.0,'nhat':50,'N':2},
    #               {'eta':0.1,'nhat':50,'N':12},
    #               {'eta':1.0,'nhat':50,'N':12},
    #               {'eta':10.0,'nhat':50,'N':12},
    #               {'eta':1.0,'nhat':1,'N':2},
    #               {'eta':1.0,'nhat':10,'N':2},
    #               {'eta':1.0,'nhat':1000,'N':2}]

    # conditionlist = [{'eta':10.0,'nhat':50,'N':2,'burst':-1,'mode':'lineage'},
    #                {'eta':1.0,'nhat':50,'N':2,'burst':-1,'mode':'lineage'},
    #                {'eta':10.0,'nhat':50,'N':12,'burst':-1,'mode':'lineage'},
    #                {'eta':1.0,'nhat':50,'N':12,'burst':-1,'mode':'lineage'}]

    conditionlist = [{'eta':1.0,'nhat':50,'N':12,'burst':-1,'mode':'lineage'}]


    # condition list used to for the line plot 
    # conditionlist = [{'nhat': 50,'N':2,'burst':-1,'mode':'lineage'},
    #                     {'nhat': 50,'N':2,'burst':-1,'mode':'population'},
    #                     {'nhat': 50,'N':2,'burst':10,'mode':'lineage'},
    #                     {'nhat': 50,'N':2,'burst':10,'mode':'population'}]

    # for condition in conditionlist:
    #   eta = condition['eta']
    #   nhat = condition['nhat']
    #   N = condition['N']

    if loaddata is False:
        if (output == 'grid'):
            f = open('grid_mode_{}_burst_{}_N.dat'.format(mode,modeburst,N),'w+') # File to write data
            gridsizex = 100
            gridsizey = 100

        elif (output == 'distribution'):
            gridsizex = 1
            gridsizey = len(conditionlist)

        elif (output == 'line'):
            gridsizex = 30
            gridsizey = len(conditionlist) # number of conditions to compare

        for inh, lognhat in enumerate(np.linspace(-1,2,gridsizey)):
            if output == 'line':
                line = np.zeros((0,2)) # Empty array with two columns to paint a line
            for logeta in np.linspace(-1,1,gridsizex):

                if output == 'grid':        
                    eta = 10**logeta
                    nhat = 10**lognhat

                if output == 'distribution':
                    eta = conditionlist[inh]['eta']

                elif output == 'line':
                    eta = 10**logeta

                if (output in ['distribution','line']):
                    nhat = conditionlist[inh]['nhat']
                    N = conditionlist[inh]['N']
                    burst = conditionlist[inh]['burst']
                    mode = conditionlist[inh]['mode']


                if burst>0:
                    modeburst = 'bursty'
                else:
                    modeburst = 'constitutive'


                W=int(N/2)

                pardict['k1'] = eta
                pardict['k0'] = nhat*eta/np.abs(burst)

                # Choosing initial condition (should not be really important)
                if burst>0:
                    pardict['m0'] = int(bt.meanpopErlang(burst,nhat*eta,N,eta,W,N,mode='bursty'))
                else:
                    pardict['m0'] = int(meantrajErlang(W/N,nhat,eta,1/N))

                pardict['mRNAgeo'] = burst
                pardict["T"] = 1

                if mode == 'lineage':
                    pardict['ffwrite'] = 1 # get the whole trajectory
                    coloridx = 0
                else:
                    pardict['ffwrite'] = 4 # get the final population
                    coloridx = 2

                # Simulations of population
                pardict['stocycle'] = 1
                pardict['presynthesis'] = W ## the actual condition is n<presynthesis
                pardict["cellphaserates"] = np.ones(N)/(1/N)
                pardict['phasenumber'] = N
                print("Single traj loop for eta={} and N ={} and nhat = {}" .format(eta,N,nhat))

                mRNA_histo_traj = np.ones(0)
                mRNA_stage_traj = np.ones(0)

                if mode == 'lineage':
                    durationfactor = 60000 # used for the grid
                    #durationfactor = 600000 # used for the line
                    correctionfactor = 10/eta # factor used to correct that distributions are noisier for lower eta
                    simulationtime = correctionfactor*durationfactor*max(1/pardict['k1'],1/pardict['k0'],1)
                    # Simulation for a trajectory of stochastic cell cycle
                    pardict['totaltime'] = simulationtime
                    pardict['ffwrite'] = 0
                    timelapse = simulationtime/(200000)
                    pardict['timelapse'] = timelapse
                    B = gill.RunTime_C(pardict)

                    simulationburnout = 0.4*simulationtime
                    sliceburnout = B[:,0]>simulationburnout

                else:

                    pardict['totaltime'] = 20
                    pardict['ffwrite'] = 4
                    B = gill.RunTime_C_Pop(pardict)
                    sliceburnout = B[:,0]>0# all the points

                bins = np.arange(min(B[sliceburnout,1]),max(B[sliceburnout,1])+2)

                skewness_simulation = skew(B[sliceburnout,1])
                kurtosis_simulation = kurtosis(B[sliceburnout,1],fisher = False) # otherwise you get excess kurtosis
                bimodality_distribution = (1+skewness_simulation*skewness_simulation)/kurtosis_simulation

                if output == 'distribution':
                    plt.hist(B[sliceburnout,1], bins= bins-0.5,width = 1.0,density = True,color = sbcolorcyclelight[1+coloridx])

                print('average from simulation', np.mean(B[sliceburnout,1]))
                print('skewness',skewness_simulation)
                print('kurtosis',kurtosis_simulation)
                print('bistable coefficient from simulation', bimodality_distribution)

                ### CREATING SUM OF NEGATIVE BINOMIALS
                # COMPUTING THE STATISTICS FOR EACH PHASE
                r_vec = nhat*eta*np.ones(N)
                k_vec = N*np.ones(N) # this is changed to get pop measurements
                if mode == 'population':
                    k_vec = k_vec*2**(1/N)
                r_vec[W:] = r_vec[W:]*2 
                d_vec = eta*np.ones(N)
                beta_vec = np.ones(N)*burst

                n0 = bt.factorial_moment(0,beta_vec,r_vec,k_vec,d_vec,mode = modeburst)
                n1 = bt.factorial_moment(1,beta_vec,r_vec,k_vec,d_vec,mode = modeburst)
                n2 = bt.factorial_moment(2,beta_vec,r_vec,k_vec,d_vec,mode = modeburst)
                n3 = bt.factorial_moment(3,beta_vec,r_vec,k_vec,d_vec,mode = modeburst)
                n4 = bt.factorial_moment(4,beta_vec,r_vec,k_vec,d_vec,mode = modeburst)

                aver = n1/n0
                sigma =  n2/n0 - aver*aver + aver
                x2 = n2/n0 + aver # second moment
                x3 = n3 + 3 *(sigma + aver*aver) - 2*aver # third moment
                x4 = n4 + 6 * x3 - 11*(sigma + aver*aver) + 6*aver # fourth moment

                if mode == 'lineage':
                        pi = 1/N*np.ones(N)
                elif mode == 'population':
                        pi = (2**(1/N)-1)/(2**(((np.arange(1,N+1)+1)/N)-1))
                aver_total = np.sum(aver*pi)
                x2_total = np.sum(x2*pi) # second moment of the whole distribution
                sigma_total = x2_total - aver_total*aver_total
                x3_total = np.sum(x3*pi)
                x4_total = np.sum(x4*pi)

                
                # negative binomial
                nb_p = aver/sigma # Note that Wikipedia and Scipy definition of the NB require  p -> 1-p 
                nb_r = aver*nb_p/(1-nb_p)

                # binomial
                # b_p = 1-sigma/aver
                # b_n = aver/b_p  
                # print('binomial parameters', b_p,b_n)

                # ADDING THE N DISTRIBUTIONS
                final_pmf_nb = np.zeros_like(bins[:-1]) # negative binomial
                final_pmf_poisson = np.zeros_like(bins[:-1]) # negative binomial
                # final_pmf_b = np.zeros_like(bins) # binomial
                for phase in range(N):
                    rv_nb = nbinom.pmf(bins[:-1],nb_r[phase],nb_p[phase])
                    rv_poisson = poisson.pmf(bins[:-1],aver[phase])
                    if mode == 'lineage':
                        pi = 1/N
                    elif mode == 'population':
                        pi = (2**(1/N)-1)/(2**(((phase+1)/N)-1))
                    # rv_b = binom.pmf(bins,b_n[phase],b_p[phase])
                    if (output == 'distribution' and showcomponents):
                        plt.plot(bins[:-1],rv_nb,'-',color = sbcolorcyclelight[2+coloridx],zorder = 1, lw = 1)


                    final_pmf_nb = final_pmf_nb + rv_nb*pi
                    final_pmf_poisson = final_pmf_poisson + rv_poisson*pi

                    # final_pmf_b = final_pmf_b + rv_b
                # final_pmf_b /= N

                #### DISTANE BETWEEN DISTRIBUTIONS
                histo_sim,edges = np.histogram(B[sliceburnout,1],bins-0.5,density = True)
                idxs  = histo_sim>0
                DKL = np.nansum(final_pmf_nb[idxs]*np.log(final_pmf_nb[idxs]/histo_sim[idxs]))
                #Cross Entropy
                #DKL = -np.nansum(final_pmf_nb[idxs]*np.log(histo_sim[idxs]))

                print('average from approximation', np.sum(final_pmf_nb*bins[:-1]))
                print('Kullback-Leibler divergence: ',DKL)
                # print('Bimodality: ',bimodality)


                if output == 'distribution':
                    plt.plot(bins[:-1],final_pmf_nb,'-',lw = 2, color = sbcolorcycledark[0+coloridx],zorder = 1)
                    plt.plot(bins[:-1],final_pmf_poisson,':',color = sbcolorcycledark[2+coloridx], zorder = 2)

                    # plt.plot(bins,final_pmf_b,'-',color = sbcolorcycledark[2], label = 'Binomial')

                    plt.xlabel('mRNA number')
                    plt.ylabel('P')
                    # plt.legend()

                    if burst<0:
                        plt.title('$\\hat n = {},  \\eta = {}, \\ N = {}, constitutive$'.format(nhat,eta,N))
                    else:
                        plt.title('$\\hat n = {},  \\eta = {}, \\ N = {}, \\beta = {}$'.format(nhat,eta,N,burst))

                    plt.savefig('NB_comp_hat{}_eta{}_N{}_{}_{}.pdf'.format(nhat,eta,N,mode,modeburst))
                    plt.clf()

                elif output == 'grid':
                    print('writing in file', "{} {} {} {} {} {}\n".format(eta,nhat,N,burst,DKL,bimodality_distribution) )
                    f.write("{} {} {} {} {} {}\n".format(eta,nhat,N,burst,DKL,bimodality_distribution))

                elif output == 'line':
                    line = np.vstack((line,[eta,DKL]))

                print('\n')

            if output == 'line':
                if mode == 'lineage':
                    color = sbcolorcyclelight[0]
                else:
                    color = sbcolorcyclelight[1]
                if burst < 0:
                    linetype = ':'
                elif burst >0:
                    linetype = '-'    

                plt.plot(line[:,0],line[:,1],linetype,color = color,label = '{} {}'.format(mode,burst))
                plt.legend()
                plt.xscale('log')
                plt.yscale('log')
                f = open('DKL_line_mode_{}_burst_{}_N_{}.dat'.format(mode,modeburst,N),'w+') # File to write data
                np.savetxt(f,line)
                f.close()


        if output == 'line':
            ax = plt.gca()
            ax.set_aspect(aspect = 0.3)
            plt.savefig('DKL_line_nhat_{}_N_{}.pdf'.format(nhat,N))


        if output == 'grid':
            f.close()

    if output == 'grid':
        A = np.loadtxt('grid_mode_{}_burst_{}_N{}.dat'.format(mode,modeburst,N ))
        x_list = np.log10(A[:,0])
        y_list = np.log10(A[:,1])

        # Computing bistable line
        x_set = np.sort(np.array(list(set(A[:,0]))))
        x_array = np.zeros(0)
        y_array = np.zeros(0)
      
        for ix,x in enumerate(x_set):
            subA = A[A[:,0]==x]
            minmax = (subA[:-1,5]-0.55)*(subA[1:,5]-0.55)
            try:
                idxs = np.where(minmax<0)[0][-1] # position of the last index
                if subA[idxs,1]>10: # to remove cases with low values of mRNA
                    y_array = np.append(y_array,subA[idxs,1])
                    x_array = np.append(x_array,x)
            except:
                pass

        dy = (max(y_list)-min(y_list))/(np.sqrt(len(A)))
        dx = (max(x_list)-min(x_list))/(np.sqrt(len(A)))

        y, x = np.mgrid[slice(min(y_list), max(y_list) + dy, dy),
                slice(min(x_list), max(x_list) + dx, dx)]

        A[:,4][A[:,4]<10**-4] = 10**-4 
        z = np.log10(A[:,4]).reshape((int(np.sqrt(len(A))),int(np.sqrt(len(A)))))
        im = plt.pcolormesh(x, y, z,vmin = -4, vmax = -1.5)
        clb = plt.colorbar(im)
        clb.set_label('$D_{KL}$', labelpad=-40, y=-0.05, rotation=0)
        plt.xlabel('$log_{10} \\eta$')
        plt.ylabel('$log_{10} \\hat{n}$')


        spl = UnivariateSpline(np.log10(x_array), np.log10(y_array))
        spl.set_smoothing_factor(0.5)
        print (x_array,spl(x_array))
        plt.plot(np.log10(x_array),spl(np.log10(x_array)),'-r')
        # plt.plot(np.log10(x_array),np.log10(y_array),'o')
        plt.savefig('grid_mode_{}_burst_{}_N{}.pdf'.format(mode,modeburst,N ))
        plt.show()


#################################################################################

def getbeforeafter(A):

    global mRNAbeforeDivide, mRNAafterDivide, protbeforeDivide, protafterDivide
    mRNAbeforeDivide = []
    mRNAafterDivide = []
    protbeforeDivide = []
    protafterDivide = []

     
    for irow,row in enumerate(A): 
        if (row[4]==1 and irow<(len(A)-3)):
            mRNAbeforeDivide.append(A[irow-1,1]) 
            mRNAafterDivide.append(A[irow,1])
            protbeforeDivide.append(A[irow-1,2]) 
            protafterDivide.append(A[irow,2])       


def averagefromtraj(A):
# This function computees the averages of mRNA and protein when the output contains all the reactions
# note that the average is not as simple as the average of entries since they are pondered by
# the duration of each state

    durations = A[1:,0]-A[:-1,0]
    ponderededmRNA = A[:-1,1]*durations
    ponderededprot = A[:-1,2]*durations
    avermRNA = sum(ponderededmRNA)/sum(durations)
    averprot = sum(ponderededprot)/sum(durations)

    return avermRNA,averprot

def averagefromtraj_timelapse(A):
# This function computees the averages of mRNA and protein when the output contains equidistant observations

    avermRNA = np.mean(A[:,1])
    averprot = np.mean(A[:,2])

    return avermRNA,averprot


def averagephasefromtraj(A):
# same as before but the output is for each phase   

    timeatphase = np.zeros(pardict['phasenumber'])
    mRNAperphase = np.zeros(pardict['phasenumber'])
    protperphase = np.zeros(pardict['phasenumber'])

    for irow,row in enumerate(A):
        if irow>0:
            cellstate = int(A[irow-1,3])
            duration = A[irow,0] - A[irow-1,0]
            timeatphase[cellstate] += duration
            mRNAperphase[cellstate] += A[irow-1,1]*duration
            protperphase[cellstate] += A[irow-1,2]*duration
    print ('mRNAperphase',mRNAperphase)
    print ('timeatphase',timeatphase)
    return timeatphase/sum(timeatphase), mRNAperphase/timeatphase, protperphase/timeatphase


def distphasefromtraj(A):
# similar to previous ones. This one returns Nphase histograms, one for each phase, each one describing
# the probability of finding i mRNAs

    maxmRNA = int(max(A[:,1]))
    maxstate = int(max(A[:,3]))
    arrayofhistograms = np.zeros((maxstate+1,maxmRNA+1)) # +1s to include zero

    #each entry will contain the time spent in that state (number of mRNAs and cell state)
    # then each histogram is obtained by normalizing each row

    for irow,row in enumerate(A):
        if irow>0:
            cellstate = int(A[irow-1,3])
            mRNAs = int(A[irow-1,1])
            duration = A[irow,0] - A[irow-1,0]
            arrayofhistograms[cellstate,mRNAs] += duration

    #arrayofhistograms = arrayofhistograms/arrayofhistograms.sum(axis=1)[:,None] # normalization of histograms
    arrayofhistograms /= sum(sum(arrayofhistograms))

    return arrayofhistograms  


def varfromtraj(A):
# This function computees the averages of mRNA and protein when the output contains all the reactions
# note that the average is not as simple as the average of entries since they are pondered by
# the duration of each state

    durations = A[1:,0]-A[:-1,0]
    ponderededmRNA = A[:-1,1]*A[:-1,1]*durations
    ponderededprot = A[:-1,2]*A[:-1,2]*durations
    mRNA2 = sum(ponderededmRNA)/sum(durations)
    prot2 = sum(ponderededprot)/sum(durations)

    avermRNA,averProt = averagefromtraj(A)

    varmRNA = mRNA2 - avermRNA*avermRNA
    varprot = prot2 - averProt*averProt

    return varmRNA,varprot


def hypoexponentialmean(k = pardict['cellphaserates']):

    print ('k',k)
    G = np.zeros(len(k))
    G[0] = 1 # this will be corrected later by normalization
    for iel in range(1,len(k)):
        G[iel] = k[iel-1]/k[iel]*G[iel-1]
    G = G/sum(G) # normalization
    print('G',G)

    r = pardict['k0']*np.ones_like(G)
    for iel in range(1,len(k)):
        if iel>=pardict['presynthesis']:
            r[iel] *= 2
    print('r',r)

    d = np.ones_like(r)*pardict['k1']

    # g[0] does not carry any meaning   
    g = G*r/(k+d)
    print ('g',g)

    f = np.zeros_like(g) # f(len(k)) does not carry any meaning
    for iel in range(0,len(k)-1):
        f[iel] = k[iel]/(k[iel+1]+d[iel+1]) 

    delta = np.ones_like(d)
    delta[0] = 1 # this value is not used in the calculations, but this is useful to define the recurrence relation
    for iel in range(1,len(k)):
        delta[iel] = delta[iel-1]*f[iel-1]
    print('delta',delta)


    theta = np.zeros_like(d) # again, theta[0] carries no meaning, setting it to zero is used to define the recurrence
    for iel in range(1,len(k)):
        theta[iel] = theta[iel-1]+g[iel]/delta[iel]
    theta = theta * delta
    print('theta',theta)

    avern1 = (2*G[0]*r[0] + k[-1]*theta[-1])/(2*(d[0]+k[0]) - k[-1]*delta[-1])
    avern = delta * avern1 + theta
    avern[0] = avern1
    avern = avern/G
    print('avern',avern)

    return avern

