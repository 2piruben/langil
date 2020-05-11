import numpy as np
import seaborn as sns
from scipy.stats import poisson,nbinom,binom,skew,kurtosis
import matplotlib.pyplot as plt
import gillespie_cyclo as gill
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


pardict = {
"Omega": 1,
"k0" : 10,  "k1" : 1, "k2" : 0, "k3" : 0, "T" : 1,
"m0" : 0, "p0" : 0, "ffwrite" : 1, "timelapse": 0.01,
"totaltime": 1000000, "dt" : 0.01, "runtype" : 0, "SEED": -1,
"runtimes": 1, "stocycle":1, "phasenumber":2, "presynthesis":1,
"cellphaserates" :[1.0] 
}

#Figure 1:
# Average number of mRNA along a single cell trajectory as a function of the cell cycle variability


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
# As derived by Casper from alternative route and confirmed by Mathematica  
    factor = 2**(1-w)*Delta*eta*nhat
    factor /= (-1+(2**Delta)+Delta*eta)
    return factor

def meanstarpopErlang(w,nhat,eta,Delta):
    factor = 2**(1-w)*nhat*eta
    factor /= (eta+np.log(2))
    return factor

# not working?
# def meanstarpopErlang(w,nhat,eta,Delta):
#   factor = 2*np.exp(eta*w)/(2*np.exp(eta)-1)
#   factor = 1+np.exp(-eta)-factor
#   factor = factor*2*nhat/(1+eta/np.log(2))
#   factor = nhat*(2-w)-factor
#   return factor


# def geo_series(r,n0,nf):
#   #sum_{j=n0}^nf r^j
#   if n0<nf:
#       return (r**n0-r**(nf+1))/(1-r)
#   else:
#       return 0


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

def varErlang(nhat,eta,W,N,mode='traj'):

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
        sum1+= pi[m-1]*factmom_Erlang(m,2,nhat,kd[m-1],W,N)/factmom_Erlang(m,0,nhat,kd[m-1],W,N)

    mE = meanErlang(nhat,eta,W,N,mode=mode)
    return sum1 + mE*(1-mE)



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


# def bigRold(eta,w,Delta):
# # corresponding to the worng formula (25) 

#   Rnum = 1 - 2*np.exp(eta) + np.exp(w*eta)
#   Rden = (1 - 2*np.exp(eta))*(1+(1+Delta*eta)**((w-1)/Delta)/( (1+Delta*eta)**(-1/Delta)-2))
#   R = 1- Rnum/Rden
#   return R


# def bigRnew(eta,w,Delta):
#   avern = meantrajErlang(w,1,eta,Delta)
#   avernstar = w + (1-w)*2 - 1/eta*(1-np.exp(-eta*(1-w))/(2-np.exp(-eta)))
#   return (avern-avernstar)/avern

def getFig1_2_3(panel='a',figure=1,loaddata=False):
    # Figs 1_2 correspond with trajectory measurements
    # Figs 2_4 correspond with population measurements

    if panel=='a':
        conditionlist = [{'N':2,'W':1,'nhat':1,'label':'$\\hat{n}=1$','mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':500,'label':'$\\hat{n}=500$', 'mRNAgeo':-1},
                         {'N':2,'W':1,'nhat':500,'label':'$\\hat{n}=500$', 'mRNAgeo':10}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)


    if panel=='b':
        conditionlist = [{'N':4,'W':1,'nhat':50,'label':'$\\Delta=1/4\\quad w=1/4$'},{'N':4,'W':2,'nhat':50,'label':'$\\Delta=1/4\\quad w=1/2$'},{'N':4,'W':3,'nhat':50,'label':'$\\Delta=1/4\\quad w=3/4$'}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)

    if panel=='c':
        conditionlist = [{'N':8,'W':2,'nhat':50,'label':'$\\Delta=1/8\\quad w=1/4$'},{'N':8,'W':4,'nhat':50,'label':'$\\Delta=1/8\\quad w=1/2$'},{'N':8,'W':6,'nhat':50,'label':'$\\Delta=1/8\\quad w=3/4$'}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)

    if panel=='d':
        conditionlist = [{'N':12,'W':3,'nhat':50,'label':'$\\Delta=1/12\\quad w=1/4$','mRNAgeo':-1},
                        {'N':12,'W':6,'nhat':50,'label':'$\\Delta=1/12\\quad w=1/2$','mRNAgeo':-1},
                        {'N':12,'W':9,'nhat':50,'label':'$\\Delta=1/12\\quad w=3/4$','mRNAgeo':-1}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)

    if panel=='c_burst': # new panel for figure 2 showing differences in error depending on burst size
        conditionlist = [{'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':1},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':10},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$', 'mRNAgeo':100},
                         {'N':2,'W':1,'nhat':50,'label':'$\\hat{n}=50$','mRNAgeo':-1}]
        etatheo = np.logspace(-2,2,100)
        etanum = np.logspace(-1.9,1.9,10)


    simulationtime = 500
    simulationburnout =  100
    if figure in [1,2]:
        repeats = 50
    elif figure in [3]:
        repeats = 5000
    elif figure in [4]:
        repeats = 25000



    for icondition,condition in enumerate(conditionlist):

        w = condition['W']/condition['N']
        Delta = 1/condition['N']

        etanum_corrected = etanum*(1+0.15*(icondition-1))

        if (condition['mRNAgeo']== -1 and figure in [1,2]): # if we are studying constitutive expression in a trajectory
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


        elif (condition['mRNAgeo'] > 0 and figure in [1,2]): # if we are studying bursty expression in a trajectory
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

        elif (condition['mRNAgeo']== -1 and figure in [3]): # if we are studying mean constitutive expression in a trajectory
            
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

        elif (condition['mRNAgeo']== -1 and figure in [4]): # if we are studying mean constitutive expression in a trajectory

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

        elif (condition['mRNAgeo']>0 and figure in [4]): # if we are studying mean constitutive expression in a trajectory
            
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

        elif (condition['mRNAgeo']>0 and figure in [3]): # if we are studying mean constitutive expression in a trajectory
            
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
            print('vartheo',vartheo)
            print('varR',varR)

        # Plotting theoretical prediction for R
        if figure in [1,3]:
            if panel == 'a':
                plt.plot(etatheo,R,'-',color = 'k',label='theory')
                # plt.plot(etatheo,mtheo,'-',color = 'k',label='theory')
                # plt.plot(etatheo,mtheo_old,':',color = 'g',label='theory')
            else:
                plt.plot(etatheo,meanR,'-',color = sbcolorcyclelight[icondition],label='theory')
                # plt.plot(etatheo,mtheo,'-',color = sbcolorcyclelight[icondition],label='theory')
        elif figure in [2,4]:
            plt.plot(etatheo,varR,'-',color = sbcolorcyclelight[icondition], label ='theory')
            # plt.plot(etatheo,vartheo,'-',color = sbcolorcyclelight[icondition], label ='theory')

        savelist = []

        # Plotting data for constitutive expression
        if loaddata:
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
                        if figure in [1,2]:
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
                        elif figure in [3,4]:
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
                        if figure in [1,2]:
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
                        elif figure in [3,4]:
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

                    # meanR = (meanmRNA-meanmRNAstar)/meanmRNA
                    # varR = (varmRNA-varmRNAstar)/varmRNA

                    # meanmRNA_error = np.std(marray)/np.sqrt(len(marray))
                    # mean_var_mRNA = np.var(marray) # variance of the mean
                    # varmRNA = np.mean(mvararray) # variance
                    # varmRNA_error = np.std(mvararray)/np.sqrt(len(marray))
                    # meanR = (np.mean(marray)-np.mean(mstararray))/np.mean(marray)
                    # varR = (np.mean(mvararray)-np.mean(mstarvararray))/np.mean(mvararray)
                    # errR =  (1 - meanR)*(np.std(marray)/np.mean(marray)+np.std(mstararray)/np.mean(mstararray))/np.sqrt(len(marray))
                    # errvarR = (1 - varR)*(np.std(mvararray)/np.mean(mvararray)+np.std(mstarvararray)/np.mean(mstarvararray))/np.sqrt(len(marray))
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

        if figure in [1,3]:

            # used fro 
            # if icondition == 0:
            #     rep = 500
            #     etanum_corrected *= 0.85
            # else:
            #     rep = 50
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


        elif figure in [2,4]:

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

            
            # plt.plot(etanum_corrected,varmRNA,'o',color = (0,0,0,0),
            #   markeredgewidth = 1.5,markersize = 4, markeredgecolor = sbcolorcyclelight[icondition],label=condition['label'])


        # Plotting data for bursty expression:
        #still in progress


        plt.legend()
        #plt.text(50,meanR[0]-icondition*0.02,condition['label'])
    
    plt.xlabel('$\\eta$')
    if figure in [2,4]:
        plt.ylabel('$R_\\sigma$')
        #plt.ylabel('$\\sigma^2$')
    else:
        plt.ylabel('$R$')
    plt.xscale('log')
    #plt.yticks([])
    ax = plt.gca()
    #ax.yaxis.set_major_formatter(ScalarFormatter())
    #ax.yaxis.set_major_formatter(NullFormatter())
    # ylin = np.linspace(0.023,0.038)
    # wopt = 2-np.sqrt(2)
    # plt.plot(wopt*np.ones_like(ylin),ylin,':',color=sbcolorcycleblue[0])
    # ax.set_yticks([0.025,0.03,0.035],[0.025,0.03,0.035])
    # plt.ylim([0.023,0.038])
    plt.tight_layout()
    plt.savefig('fig{}{}.pdf'.format(figure,panel))
    plt.show()

def getFig_sidehisto(output = 'trajhistogram'):

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


def getFig_trajcomparison():

    # Figure 1
    # fig = plt.figure(1, figsize=(15,8))
    # gs = gridspec.GridSpec(3,1)

    # Figure 2
    fig = plt.figure(1, figsize=(15,10))
    gs = gridspec.GridSpec(2,1)


    gs.update(left=0.05, right=0.95, bottom=0.08, top=0.93, wspace=0.02, hspace=0.00)

    nhat = 50
    # pardict['k0'] = 10.0
    tau = 1 # total cell cycle duration
    pardict["T"] = tau

    # Figure 2
    Ns = [2,12]
    etas = [0.1,5]
    Ws = [1,3]
    lims = [[0,30],[0,150]]

    # Figure 1
    # Ns = [4,4,4]
    # Ws = [3,3,3]
    # etas = [0.1,1,10]
    # lims = [[0,20],[0,100],[0,150]]

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

    # Labels Figure 1
    #ax_trajs[2].set_xticks([0.6666,1,2,3,4,5]) 
    #ax_trajs[2].set_xticklabels(['W','T','2T','3T','4T','5T'])

    # Labels Figure 2
    ax_trajs[1].set_xticks([2,4,6,8]) 
    ax_trajs[1].set_xticklabels(['2T','4T','6T','8T'])

    plt.xlabel('time')
    ax_trajs[1].set_ylabel('mRNA number')

    plt.savefig('trajseta_Fig2.pdf')
    plt.show()


def getFigDistributionComparison(mode = 'lineage', burst = -1, output = 'grid', loaddata = False, N = 12, showcomponents = False):

    # mode = lineage / population
    # output = grid / line / distribution



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



def getSchwancomparison():

    Sdf = pd.read_csv('Schwanhausser.csv')

    paper_mRNA = 'mRNA copy number experiment [molecules/cell]'
    mRNA_half_life = 'mRNA half-life experiment [h]'
    paper_k0 = 'transcription rate (vsr) experiment [molecules/(cell*h)]'

    T = 27.5
    Delta = 1/4
    degradation = np.log(2)/Sdf[mRNA_half_life]
    eta = degradation*T
    eta_cont = np.logspace(np.log10(0.4),np.log10(20))

    w = 0.5
    
    Schwan_k0 = Sdf[paper_mRNA]*degradation/(1+1/(eta)*(1-(1-np.exp(-eta))/(2-np.exp(-eta)))*(np.exp(-eta)-1))
    Schwan_k0_cont = eta_cont/T/(1+1/(eta_cont)*(1-(1-np.exp(-eta_cont))/(2-np.exp(-eta_cont)))*(np.exp(-eta_cont)-1)) # same but for a cont itnuous eta set

    theo_mRNA_lin = meantrajErlang(w,1,eta,Delta) # paper plus noise in cell cycle
    theo_mRNA_pop = meanpopErlang(w,1,eta,Delta) # paper plus noise in cell cycle

    #####################
    w = 1.0
    theo_mRNAstar_lin = meantrajErlang(w,1,eta,Delta) # this is the same they do in the paper
    theo_mRNAstar_lin_cont = meantrajErlang(w,1,eta_cont,Delta) # same but with a continuous of eta
    r_theo =  Sdf[paper_mRNA]*degradation/theo_mRNAstar_lin*(2-w)
    r_theo_cont =  eta_cont/T/theo_mRNAstar_lin_cont*(2-w)
    
    error = (Schwan_k0 - r_theo)/(Schwan_k0)
    error_cont = (Schwan_k0_cont - r_theo_cont)/(Schwan_k0_cont)
    error1 = error
    error1_cont = error_cont 
    #error1 = error[~np.isnan(error)]

    #####################
    w = 1.0
    theo_mRNAstar_pop1 = meanpopErlang(w,1,eta,Delta) # this is the same they do in the paper
    theo_mRNAstar_pop1_cont = meanpopErlang(w,1,eta_cont,Delta) # this is the same they do in the paper
    r_theo =  Sdf[paper_mRNA]*degradation/theo_mRNAstar_pop1*(2-w)
    r_theo_cont =  eta_cont/T/theo_mRNAstar_pop1_cont*(2-w)
    error = (Schwan_k0 - r_theo)/(Schwan_k0)
    error_cont = (Schwan_k0_cont - r_theo_cont)/(Schwan_k0_cont)
    error2 = error
    error2_cont = error_cont
    #error2 = error[~np.isnan(error)]

    #####################
    w = 0.5
    theo_mRNAstar_pop2 = meanpopErlang(w,1,eta,Delta) # this is the same they do in the paper
    theo_mRNAstar_pop2_cont = meanpopErlang(w,1,eta_cont,Delta) # this is the same they do in the paper
    r_theo =  Sdf[paper_mRNA]*degradation/theo_mRNAstar_pop2*(2-w)
    r_theo_cont =  eta_cont/T/theo_mRNAstar_pop2_cont*(2-w)
    error = (Schwan_k0 - r_theo)/(Schwan_k0)  
    error_cont = (Schwan_k0_cont - r_theo_cont)/(Schwan_k0_cont)
    error3 = error
    error3_cont = error_cont
    #error3 = error[~np.isnan(error)]


###### VIOLIN PLOT
    # print ('error1',error1)
    # errors = [error1,error2,error3]
    # vplot = plt.violinplot(errors, [0,1,2],
    #     showmeans=False, showmedians=False,
    #     showextrema=False)
    # for iel,el in enumerate(vplot['bodies']):
    #     el.set_facecolor(sbcolorcyclelight[iel])
    #     el.set_edgecolor(sbcolorcycledark[iel])
    #     el.set_alpha(1)
    #     q1,q2,q3 = np.percentile (errors[iel],[25,50,75])
    #     plt.scatter([iel], q2, marker='o', color='white', s=30, zorder=3)
    #     plt.vlines([iel],[q1], [q3], linestyle='-',
    #         lw=5,zorder = 2, color = sbcolorcycledark[iel])

###### HISTOGRAM PLOT (SAME AS VIOLIN)
    # print ('error1',error1)
    errors = [error1,error2,error3]
    errors_cont = [error1_cont,error2_cont,error3_cont]

    for ierror,error in enumerate(errors):
        color = list(sbcolorcyclelight[ierror])
        color.append(0.6) # alpha channel
        hplot = plt.hist(error,
            color = color,
            edgecolor = sbcolorcycledark[ierror])
    plt.ylabel('Number of genes')
    plt.xlabel('Transcription rate error')
    plt.ylim([0,1500])
    plt.xlim([-0.17,0.05])
    plt.savefig('Schwan_histo.pdf')
    plt.show()
    plt.clf()

    #     showmeans=False, showmedians=False,
    #     showextrema=False)
    # for iel,el in enumerate(vplot['bodies']):
    #     el.set_facecolor(sbcolorcyclelight[iel])
    #     el.set_edgecolor(sbcolorcycledark[iel])
    #     el.set_alpha(1)
    #     q1,q2,q3 = np.percentile (errors[iel],[25,50,75])
    #     plt.scatter([iel], q2, marker='o', color='white', s=30, zorder=3)
    #     plt.vlines([iel],[q1], [q3], linestyle='-',
    #         lw=5,zorder = 2, color = sbcolorcycledark[iel])



#### SCATTER PLOT
    # fig, ax = plt.subplots()
    # for ierror,error in enumerate(errors):
    #     color_marker = list(sbcolorcyclelight[ierror])
    #     color_edge = list(sbcolorcycledark[ierror])
    #     color_marker.append(0.1) # alpha channel
    #     color_edge.append(0.2) # alpha channel
    #     ax.scatter(eta, error,
    #               c = [color_marker], edgecolors = [color_edge],
    #               s = 20)

    #     ax.plot(eta_cont, errors_cont[ierror],'k:')

    #     # plt.scatter(eta, error2,
    #     #         c = sbcolorcyclelight[0], edgecolors = sbcolorcycledark[0], alpha = 0.2, s = 10)

    # # #xline = np.linspace(np.min(eta)*0.9,np.max(eta)*1.1)
    # # xline = np.linspace(np.min(Schwan_k0)*0.9,np.max(Schwan_k0)*1.1)
    # # yline = np.logspace(-1.5, 3)
    #     #plt.plot(xline,np.zeros_like(xline),'k:')
    # #plt.yscale('log')
    # fontsize = 16
    # plt.xscale('log')
    # ax.set_xlabel('$\\eta$', fontsize = fontsize)
    # ax.set_ylabel('Transcription rate error', fontsize = fontsize)
    #     #plt.title('population (w=0.5) vs.reproduced', fontsize = fontsize)

    # plt.ylim([-0.17,0.05])
    # #plt.xlim([eta_cont[0],eta_cont[-1]])
    # plt.tight_layout()
    # plt.savefig('Schwan_scatter_eta.pdf')
    # plt.show()

def getFig4(panel = 'a',loaddata = True):

    if panel == 'a':
        Nlist = [2,8,12]
        for i,N in enumerate(Nlist):
            pardict['phasenumber'] = N
            pardict['presynthesis'] = N/2
            pardict['stocycle'] = 1
            pardict['k1'] = eta
            pardict['k0'] = eta*nhat
            pardict['m0'] = nhat*(1+w)
            tau = 1 # total cell cycle duration
            pardict["T"] = tau
            pardict["cellphaserates"] = np.ones(N)/(tau/N)
            pardict['totaltime'] = 3
            repeats = 500
            pardict['ffwrite'] = 4
            sumprob = np.zeros(N)
            sumprob2 = np.zeros(N) 

            for repeat in tqdm(range(repeats)):
                B = gill.RunTime_C_Pop(pardict)
                prob = np.histogram(B[:,3],np.arange(N+1)-0.5)[0]
                prob = prob/np.sum(prob) # Normalizations
                sumprob += prob
                sumprob2 += prob*prob

            averprob = sumprob/repeats
            print('errorprob',sumprob2)
            averprob2 =sumprob2/repeats-averprob*averprob
            print('errorprob',averprob2)
            errorprob = np.sqrt(averprob2)/np.sqrt(repeats) 
            print('errorprob',errorprob)            

            phases = np.arange(N)+1
            #plt.scatter(phases,(2**(1/N)-1)/(2**(phases/N-1)),marker='s', s = 40, c = [(0,0,0,0)], edgecolor=sbcolorcycledark[i])
            plt.step(phases-0.5,(2**(1/N)-1)/(2**(phases/N-1)),color=sbcolorcyclelight[i],where = 'post')
            plt.step(phases+0.5,(2**(1/N)-1)/(2**(phases/N-1)),color=sbcolorcyclelight[i],where = 'pre')
            plt.scatter(phases,averprob,marker ='o',s = 10,color = sbcolorcycledark[i])
            plt.errorbar(phases,averprob,1.5*errorprob,ls='none', color =sbcolorcycledark[i],zorder = 100 )

            cont_time = np.linspace(0,N,100)
            pi_cont = np.log(2)/N*2**(1-cont_time/N)
            cont_time_adapted = cont_time+0.5
            plt.plot(cont_time_adapted,pi_cont,':',color = sbcolorcycledark[i])
        plt.xlim([0.5,12.5])
        plt.ylim([0.05,0.4])
        plt.yscale('log')
        plt.xlabel('cell phase')
        plt.ylabel('probability')
        plt.tight_layout()
        plt.savefig('fig4panel{}.pdf'.format(panel))
        plt.show()

    if (panel == 'b'):

        #etanum = np.logspace(-2,2,5)
        etanum = np.logspace(-2,2,6)
        # etanum = np.array([2.0])
        etatheo = np.logspace(-2,2,100)
        Nlist = [2]
        nhat = 50
        W = 1

        for icondition,N in enumerate(Nlist):
            if loaddata:
                D = np.loadtxt('fig4_panelb_cond{}.dat'.format(icondition))
            else:
                etanum_corrected = etanum*(1+0.15*(icondition-1))

                savelist = []
                for eta in etanum_corrected:
     
                    pardict['phasenumber'] = N
                    pardict['presynthesis'] = W ## the actual condition is n<presynthesis
                    pardict['stocycle'] = 1
                    pardict['k1'] = eta
                    # pardict['k1'] = 0.000001 # because cell cucle lenght is T =1
                    pardict['k0'] = nhat*eta
                    # pardict['k0'] = 10.0
                    pardict['m0'] = int(2*meantrajErlang(0.5,nhat,eta,1/N))
                    tau = 1 # total cell cycle duration
                    pardict["T"] = tau
                    pardict["cellphaserates"] = np.ones(N)/(tau/N)


                    ## Calculation of the average mRNA trajectory
                    repeats = 2
                    pardict['ffwrite'] = 0
                    #simulationtime = 600*max(4/pardict['k1'],4/pardict['k0'],tau)
                    simulationtime = 600
                    simulationburnout = 0.4*simulationtime
                    pardict['totaltime'] = simulationtime
                    timelapse = simulationtime/10000
                    pardict['timelapse'] = timelapse
                    marray = np.ones(0)         
                    mvararray = np.ones(0)          
                    print("Single traj loop for eta={} and N ={}".format(eta,N))
                    for repeat in tqdm(range(repeats)):

                        B = gill.RunTime_C(pardict)
                        sliceburnout = B[:,0]>simulationburnout
                        m,p = averagefromtraj_timelapse(B[sliceburnout])
                        mvar,pvar = varfromtraj(B[sliceburnout])
                        marray = np.append(marray,m)
                        mvararray = np.append(mvararray,mvar)
                    meanmarray = np.mean(marray)
                    varmarray = np.mean(mvararray)
                    errormarray = np.std(marray)/len(marray) 
                    errorvarmarray = np.std(mvararray)/len(marray) 


                    # Calculation of the average mRNA population
                    marray_pop = np.ones(0)         
                    mvararray_pop = np.ones(0)  
                    pardict['ffwrite'] = 4
                    pardict['totaltime'] = 12
                    print("Pop traj loop for eta={} and N={}".format(eta,N))
                    repeats = 5
                    for repeat in tqdm(range(repeats)):
                        B = gill.RunTime_C_Pop(pardict)
                        marray_pop = np.append(marray_pop,np.mean(B[:,1]))
                        mvararray_pop = np.append(marray_pop,np.var(B[:,1]))
                    meanmarray_pop = np.mean(marray_pop)
                    errormarray_pop = np.std(marray_pop)/len(marray_pop) 
                    varmarray_pop = np.mean(mvararray_pop)
                    errorvarmarray_pop = np.std(mvararray_pop)/len(marray_pop) 

                    R = (meanmarray-meanmarray_pop)/(meanmarray)
                    errR =  (1 - R)*(errormarray+errormarray_pop)/np.sqrt(len(marray))
                    varR = (varmarray-varmarray_pop)/(varmarray)
                    errvarR =  (1 - varR)*(errorvarmarray+errorvarmarray_pop)/np.sqrt(len(marray))
                    savelist.append([eta,meanmarray,errormarray,meanmarray_pop,errormarray_pop,R,errR,varR,errvarR])
                D = np.array(savelist)  
                np.savetxt('fig4_panelb_cond{}.dat'.format(icondition),D)
            
            #plt.scatter(D[:,0],D[:,5]) ## R
            #plt.scatter(D[:,0],D[:,1]) ## meantraj
            plt.scatter(D[:,0],D[:,3]) ## meanpop


            Rtheo = (meantrajErlang(W/N,nhat,etatheo,1/N) - meanpopErlang(W/N,nhat,etatheo,1/N))/meantrajErlang(W/N,nhat,etatheo,1/N)
            mpoptheo = meanpopErlang(W/N,nhat,etatheo,1/N)
            mtrajtheo = meantrajErlang(W/N,nhat,etatheo,1/N)

            plt.plot(etatheo,mpoptheo,'-',color = sbcolorcycledark[icondition])
            plt.xscale('log')
        # Plotting deterministic case
        #Rtheo = (meantrajErlang(0.5,nhat,etatheo,1E-8) - meanpopErlang(0.5,nhat,etatheo,1E-8))/meantrajErlang(0.5,nhat,etatheo,1E-8)
        #plt.plot(etatheo,Rtheo)
        plt.xlabel('$\\eta$')
        plt.ylabel('$(n^p)_1$')
        plt.savefig('fig4_panel{}.pdf'.format(panel))
        plt.show()





def getFig_testburst(panel = 'a',loaddata = False):

    if panel == 'a':

        conditionlist = [{'N':2,'W':1,'nhat':1,'beta':5,'label':'$n = 1$'},
        {'N':2,'W':1,'nhat':50,'beta':5,'label':'$n = 50$'},
        {'N':2,'W':1,'nhat':500,'beta':5,'label':'$n = 500$'},
        ]

    if panel == 'b':

        conditionlist = [{'N':12,'W':3,'nhat':50,'beta':5,'label':'$w = 1/4$'},
        {'N':12,'W':6,'nhat':50,'beta':5,'label':'$w = 1/2$'},
        {'N':12,'W':9,'nhat':50,'beta':5,'label':'$w = 3/4$'},
        ]

    if panel == 'c':

        conditionlist = [{'N':2,'W':1,'nhat':50,'beta':1,'label':'$\\beta = 1$'},
        {'N':2,'W':1,'nhat':50,'beta':5,'label':'$\\beta = 5$'},
        {'N':2,'W':1,'nhat':50,'beta':10,'label':'$\\beta = 10$'},
        ]



    etatheo = np.logspace(-2,2,100)
    etanum = np.logspace(-1.9,1.9,10)
    # etanum = np.array([1])

    for icondition,condition in enumerate(conditionlist):

        w = condition['W']/condition['N']
        nhat = condition['nhat']
        Delta = 1/condition['N']

        repeats = 200

        etanum_corrected = etanum*(1+0.15*(icondition-1))

        #mstartheo = meanstartrajErlang(w,nhat,etatheo,Delta)
        #mtheo = meantrajErlang(w,nhat,etatheo,Delta)
        mtheo = np.ones(0)
        mstartheo = np.ones(0)
        maxN = 1000 # N used to get the deterministic limit
        for eta in etatheo:
            mtheo = np.append(mtheo,bt.meanErlang(condition['beta'],nhat*eta,1/Delta,eta,condition['W'],condition['N']))
            mstartheo = np.append(mstartheo,bt.meanErlang(condition['beta'],nhat*eta,maxN,eta,int(condition['W']/condition['N']*maxN),maxN))
        R = (mtheo-mstartheo)/mtheo


        # varstartheo = varstartrajErlang(condition['nhat'],etatheo,w,Delta)
        print(condition['nhat'],etatheo,w,condition['W'],Delta)
        #vartheo = vartrajErlang(condition['nhat'],etatheo,condition['W'],Delta)
        vartheo = np.ones(0)
        varstartheo = np.ones(0)
        for eta in etatheo:
            vartheo = np.append(vartheo,bt.varErlang(condition['beta'],nhat*eta,1/Delta,eta,condition['W'],condition['N']))
            varstartheo = np.append(varstartheo,bt.varErlang(condition['beta'],nhat*eta,maxN,eta,int(condition['W']/condition['N']*maxN),maxN))
        varR = (vartheo-varstartheo)/vartheo
        # print('varstarteo',varstartheo)
        # print('vartheo',vartheo)
        # print('varR',varR)

        plt.plot(etatheo,varR,'-',color=sbcolorcyclelight[icondition])
        # plt.plot(etatheo,vartheo,'-',color = 'k',label='bursty')
        savelist = []

        # Plotting data for constitutive expression
        if loaddata:
            D = np.loadtxt('burst_panel_{}{}.dat'.format(panel,icondition))
        else:
            for eta in etanum_corrected:
                    # simulation time
                    pardict['phasenumber'] = condition['N']
                    pardict['presynthesis'] = condition['W']
                    pardict['stocycle'] = 1
                    pardict['mRNAgeo'] = condition['beta']
                    pardict['k1'] = eta # because cell cucle lenght is T =1
                    if pardict['mRNAgeo']>0:
                        pardict['k0'] = condition['nhat']*eta/pardict['mRNAgeo']
                    else:
                        pardict['k0'] = condition['nhat']*eta
                    # pardict['m0'] = condition['nhat']
                    pardict['m0'] = 0
                    tau = 1 # total cell cycle duration
                    pardict["T"] = tau
                    pardict["cellphaserates"] = np.ones(condition['N'])/(tau/condition['N'])


                    simulationtime = 600*max(4/pardict['k1'],4/pardict['k0'],tau)
                    timelapse = simulationtime/100000
                    # simulationtime = 10
                    # timelapse = 0.1
                    simulationburnout = 0.4*simulationtime
                    pardict['totaltime'] = simulationtime
                    pardict['timelapse'] = timelapse
                    pardict['ffwrite'] = 0
                    print("Calculating condition eta:{}, condition:{}, simulationtime:{}, timelapse:{}".format(eta,icondition,simulationtime,timelapse))

                    marray = np.ones(0)         
                    mvararray = np.ones(0)          

                    for repeat in tqdm(range(repeats)):

                        B = gill.RunTime_C(pardict)
                        sliceburnout = B[:,0]>simulationburnout
                        m,p = averagefromtraj_timelapse(B[sliceburnout])
                        mvar,pvar = varfromtraj(B[sliceburnout])
                        marray = np.append(marray,m)
                        mvararray = np.append(mvararray,mvar)

                    # Block to compute the deterministic computationally
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
                        B = gill.RunTime_C(pardict)
                        sliceburnout = B[:,0]>simulationburnout
                        m,p = averagefromtraj_timelapse(B[sliceburnout])
                        mvar,pvar = varfromtraj(B[sliceburnout])
                        mstararray = np.append(mstararray,m)
                        mstarvararray = np.append(mstarvararray,mvar)

                    meanR = (np.mean(marray)-np.mean(mstararray))/np.mean(marray)
                    varR = (np.mean(mvararray)-np.mean(mstarvararray))/np.mean(mvararray)
                    errR =  (1 - meanR)*(np.std(marray)/np.mean(marray)+np.std(mstararray)/np.mean(mstararray))/np.sqrt(len(marray))
                    errvarR = (1 - varR)*(np.std(mvararray)/np.mean(mvararray)+np.std(mstarvararray)/np.mean(mstarvararray))/np.sqrt(len(R))
                    savelist.append([eta,meanR,errR,varR,errvarR])  
                    # savelist.append([eta,np.mean(marray),np.var(marray),np.std(marray)/np.sqrt(len(marray)),
                                    # np.mean(mvararray),np.var(mvararray),np.std(mvararray)/np.sqrt(len(mvararray))])  
            # np.savetxt('fig1_panel{}_cond{}.dat'.format(panel,icondition),np.array(savelist))s
            np.savetxt('burst_panel_{}{}.dat'.format(panel,icondition),savelist)
            D = np.array(savelist)
        plt.plot(D[:,0],D[:,3],'o',color = (0,0,0,0),
                markeredgewidth = 1.5,markersize = 4, markeredgecolor = sbcolorcyclelight[icondition],label=condition['label'])
        plt.errorbar(D[:,0],D[:,3],yerr =1.5*D[:,4],ls='none',color=sbcolorcyclelight[icondition])#fmt=None,color = sbcolorcycleblue[icondition])
        
        # Plotting data for bursty expression:
        #still in progress


        plt.legend()
        #plt.text(50,meanR[0]-icondition*0.02,condition['label'])
    
    plt.xlabel('$\\eta$')
    plt.ylabel('$\\sigma^2$')
    plt.xscale('log')
    #plt.yticks([])
    #ax.yaxis.set_major_formatter(ScalarFormatter())
    #ax.yaxis.set_major_formatter(NullFormatter())
    # ylin = np.linspace(0.023,0.038)
    # wopt = 2-np.sqrt(2)
    # plt.plot(wopt*np.ones_like(ylin),ylin,':',color=sbcolorcycleblue[0])
    # ax.set_yticks([0.025,0.03,0.035],[0.025,0.03,0.035])
    # plt.ylim([0.023,0.038])
    plt.tight_layout()
    plt.savefig('figbusrst_{}.pdf'.format(panel))
    plt.show()


def getHistoDecays():

    eta_fibroblast = np.array([0])
    decay_yeast = np.array([0])
    decay_ecoli = np.array([0])

    # load data from fibroblast
    Sdf = pd.read_csv('Schwanhausser.csv')

    paper_mRNA = 'mRNA copy number experiment [molecules/cell]'
    mRNA_half_life = 'mRNA half-life experiment [h]'
    paper_k0 = 'transcription rate (vsr) experiment [molecules/(cell*h)]'

    T = 27.5
    Delta = 1/4
    degradation = np.log(2)/Sdf[mRNA_half_life]
    eta_fibroblast = degradation*T

    # load data from ecoli
    Sdf = pd.read_csv('halflifeBacteria.txt', sep = '\t', skiprows = 1)
    degradation = np.log(2)/Sdf['hl']
    T = 20.5 # Doubling time in minutes From (Liang et al JMB 1999)
    eta_ecoli = degradation * T

    # load data from yeast
    Sdf = pd.read_csv('halflifeYeast.csv', skiprows = 1)  
    degradation = Sdf['decay constant (m)']     
    T = 90
    eta_yeast = degradation * T


    for ierror,eta in enumerate([eta_fibroblast,eta_yeast,eta_ecoli]):
        print(ierror, 'stats:', np.mean(eta),np.std(eta))
        lenbins = np.sqrt(len(eta))
        logbins = np.logspace(np.log10(0.1),np.log10(100),lenbins)

        color = list(sbcolorcyclelight[ierror])
        color.append(0.4) # alpha channel
        hplot = plt.hist(eta,
            histtype = 'stepfilled',
            bins=logbins,
            density = True,
            color = color,
            edgecolor = sbcolorcycledark[ierror])
    plt.xscale('log')
    plt.ylabel('Density')
    plt.xlabel('$\\eta$')

    #plt.ylim([0,1500])
    #plt.xlim([-0.17,0.05])
    plt.savefig('eta_species.pdf')
    plt.show()
    plt.clf()







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

def checkstatistics():

    A = np.loadtxt('output_gillespie_cyclo.out')
    
    print("Average final mRNA number is {}, expected number is {}".format(np.mean(A[:,5]),50*2.0*(1-1/(2-np.exp(-1.7))*np.exp(1.7/2)  )))
    print("pi2 compurational:{} , theoretical:{}".format(np.mean(A[:,3]),np.sqrt(2)-1 ))
    print("Average time inside the cell cycle {}, theoretical: {}".format(np.mean(A[:,4]), 2**(-3/2)))

    gamma = 2/1.7
    deltafrac = (gamma/(1+gamma))
    thetaN = 50/2*(2-2*deltafrac)
    n11 = (2*50+2*gamma*thetaN)/(2*2*(1+gamma)-2*gamma*deltafrac)
    npart = deltafrac*n11+50/2*(2-2*deltafrac)


    print("Average initial mRNA number is {}, expected number is {}".format(np.mean(A[:,5])/2,npart))


# eta, delta, r = sp.symbols('eta delta r')
# sp.init_printing(use_unicode=True)

# print(sp.limit(        1 -     (r + 2*(1-r) - 1/eta*(1-sp.exp(-eta*(1-r))/(2-sp.exp(-eta))))  / (r + 2*(1-r) - 1/eta*(1-(1+eta*delta)**((r-1)/delta)/(2-(1+eta*delta)**(-1/delta))))   ,eta,0))
