#! /bin/sh
""":"
exec python $0 ${1+"$@"}
"""
#"
###################################################################################
###################################################################################
###         SCRIPT FOR RUNNING HIV langil
################################################# Ruben Perez-Carrasco

# This module creates a set of input files readable by the HIVlangil program
# and writes the results to a certain output file. It allows to create
# several input files with different parameters to run them. The files are created in
# a folder with the name of the day the script was created and makes a backup
# of anyfile that would be overwritten by mistake when running a program with
# the same input/output name

# Additionally the script is able to create a pool of processes so different cALls
# for different input files can be run in parallel, one in each processor


from __future__ import print_function
from string import Template
from subprocess import call
import numpy as np
import os
import shutil
import time
import datetime
from multiprocessing import Process,Pool,Lock  #parallel processing
import multiprocessing as mp

pars= dict ({'Omega':5000,'lambdaT':10000, 'deltaT':0.0166,
            'k':2.4e-8,'eps':0.85,'eta':0.001,
            'aL':0.1,'aLs':0.1,'aLC':0.0,'aLC2':0.0,'aLT':1.0,
            'deltaL':0.001,'r':0.2,'Lmax':1.888,
            'deltaTS':1.0,
            'lambdaV':2000,'deltaV':23,
            'adiabaticmethod':4})

pars.update({'timelapse':1,'totaltime':100.0,'SEED':-1,'ffwrite':0,
            'dt':0.001,'runtype':3})
# runtype= {GILLESPIE=0,CLE=1,MACROSCOPIC=2, HEUN = 3}
# ffwrite=3 TIMESTEPMULT: record every time following formula T(i)=dt*(timelapse^i-1) 

# Mean First Passage Time parallelepipeds
pars.update({'Tf':599999,'Lf':0.934778,'TSf':0.115067,'Vf':10,
            'T0':600000,'L0':0.01,'TS0':0.01,'V0':1})

pars.update({'DT0':60000,'DL0':0.09,'DTS0':0.01,'DV0':1,
            'DTf':60000,'DLf':0.09,'DTSf':0.01,'DVf':1})

# Definition of the templates that will be used for the output and input files and will be
# subsituted by the dictionariy already defined

#paramstringtemplate & simstringtemplate are headers with the summary of the parameters used
paramsstringtemplate = Template('''# Parameters used:
# Omega=$Omega,'
# lambdaT=$lambdaT, deltaT=$deltaT,
# k=$k', eps=$eps, eta=$eta,
# aL=$aL, aLs=$aLs, aLC=$aLC, aLC2=$aLC2, aLT=$aLT,
# deltaL= $deltaL, r= $r, Lmax: $Lmax,
# adiabaticmethod=$adiabaticmethod,
# deltaTS=$deltaTS, lambdaV=$lambdaV, deltaV=$deltaV
# T0=$T0, L0=$L0, TS0=$TS0, V0=$V0''')

simstringtemplate = Template('''# With simulation parameters:
#   timelapse=$timelapse, totaltime=$totaltime, SEED=$SEED,ffwrite=$ffwrite. dt=$dt, runtype=$runtype ''')

inputfiletemplate = Template('''# Control Parameters

Omega $Omega # System Volume

# Parameters of the model

lambdaT $lambdaT // production and degradation of healthy cells
deltaT $deltaT
k $k // infection rate
eps $eps // drug efficacy
eta $eta // fraction of inections resulting in L

deltaL $deltaL // death rate of L
r $r // L homeostasis speed
Lmax $Lmax // homoeostatic carrying capacity

aL $aL // rate of activation from latency (will depend on parameters bellow)
aLs $aLs
aLC $aLC 
aLC2 $aLC2
aLT $aLT // parameters controlling aL in time

adiabaticmethod $adiabaticmethod // number of stochastic reactions

deltaTS $deltaTS; // death of active infected cells
lambdaV $lambdaV // virus production and clearance rate 
deltaV $deltaV 

# Initial and Final Conditions (in concentration units)

T0 $T0
L0 $L0
TS0 $TS0
V0 $V0

Tf $Tf
Lf $Lf
TSf $TSf
Vf $Vf

# pass jump transition recording

DT0 $DT0
DL0 $DL0
DTS0 $DTS0
DV0 $DV0
DTf $DTf
DLf $DLf
DTSf $DTSf
DVf $DVf

# Parameters for the integration

ffwrite $ffwrite # kind of recording protocol
#                0 TIMESTEP: record every timelapse time span
#                1 ALL: record after every reaction(integration step) takes place
#                2 NOWRITE:
#                3 TIMESTEPMULT: record every time following formula T(i)=dt*(timelapse^i-1) 
#                              note that is required timelapse>1.0

timelapse $timelapse # time between recordingds if fwrite==0
totaltime $totaltime
dt  $dt #Integration time step for Gillespie and Deterministic && is dx for CLE
runtype $runtype Run_Type {GILLESPIE=0,CLE=1,MACROSCOPIC=2}
SEED $SEED''')

# Creation of a folder with the name of the current date

day=time.strftime("%d%b%y")
outputdirname="output"+day
if not os.path.exists("output"+day): # if the folder doesn't exist create it
    os.makedirs(outputdirname)
outputdirname=outputdirname+'/'

# Creation of a log file. It will contain info of the running script along the day

logname = outputdirname+"HIVlangil"+day+".log"
with open(logname,"a") as logfile:
    logfile.write("#"*43+'\n')
    logfile.write("#"*43+'\n')
    logfile.write(time.strftime("Starting script at %c \n"))
    logfile.write("#"*43+'\n')
    logfile.write("#"*43+'\n\n')

# different cALls to different programs will be first stored in a an list of dictionaries cALl
calldict={"shellorder":"", "basename":"","loglocalname":""}
calldicts=[] # list containing the cALls to map it into the process pool

########### LOOP ALONG THE DIFFERENT PARAMETERS WANTED
batch=[
#{"runtype":2,"label":"MACRO","N":1,"adiabaticmethod":4},
#{"runtype":1,"label":"CLE","N":10,"adiabaticmethod":4},
{"runtype":3,"label":"GILL1","N":10,"adiabaticmethod":1}]
#{"runtype":0,"label":"GILL2","N":10,"adiabaticmethod":2},
#{"runtype":0,"label":"GILL3","N":10,"adiabaticmethod":3},
#{"runtype":0,"label":"GILL4","N":1,"adiabaticmethod":4}]

for exper in batch:
    pars["adiabaticmethod"]=exper["adiabaticmethod"]
    pars["runtype"]=exper["runtype"]
    for N in np.arange(exper["N"]): # deeper loop over different parameters
        basename=(exper["label"]+"N"+str(N)) # name of the files
        inputname=outputdirname+basename+".in"
        outputname=outputdirname+basename+".out"
        loglocalname = outputdirname+basename+".log"
        if os.path.exists(outputname): # if outputfile exists make a backup instead of  overwriting
            outputnameold=outputname+time.strftime("%d%b%H%M")
            shutil.copyfile(outputname,outputnameold)
        with open(loglocalname,"a") as loglocalfile:
            loglocalfile.write("Inputfile created at: "+time.strftime("%X %a %d %b '%y \n"))
            loglocalfile.write(paramsstringtemplate.substitute(pars)+'\n')
        with open(inputname,"w") as inputfile:
            inputfile.write(inputfiletemplate.substitute(pars))
                    # Shell order contains the shell command to be run
        shellorder="./HIVlangil "+inputname+" "+outputname+" >>"+loglocalname
        calldict["shellorder"]=shellorder
        calldict["basename"]=basename
        calldict["loglocalname"]=loglocalname
        calldicts.append(dict(calldict)) # dict is mandatory to append by value not by  reference
            #               p=Process(target=cALlmain, args=(shellorder,)) # parallel code
            #               p.start()

            # Definition of the script section being run in different processors "cALlmain"
        lock=Lock() # lock for writing safely into logfile

def callmain(calldict):

    lock.acquire()
    with open(logname,"a") as logfile:
        logfile.write("Starting process "+calldict["basename"]+" at "+time.strftime("%X %x")+"\n")
    lock.release()

    initialtime=time.time()
    call(calldict["shellorder"], shell=True) # cALl to one cALl inside cALls

    lock.acquire()
    with open(logname,"a") as logfile:
        logfile.write("Finished process "+calldict["basename"]+" ! at "+time.strftime("%X %a %d %b '%y")+"\n")
        duration=time.time()-initialtime
        logfile.write("Duration is "+str(duration)+" seconds ("+
                     str(datetime.timedelta(seconds=duration))+")\n")
        lock.release()

if __name__ == '__main__': # This is necessary... (don't ask me why)

    cpunum=mp.cpu_count()
    pool = Pool(processes=cpunum) # creating a pool with processors equal to the number of processors
    print("Number of processors "+str(cpunum))
    pool.map(callmain,calldicts)  # mapping of all the cALls necessary into the cALling function



#       time.strftime("Finishing script at %c")+"\n\n")
