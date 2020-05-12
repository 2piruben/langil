import os
import numpy as np
import matplotlib.pyplot as plt 

plt.ioff()

basedictionary = {} # input_dic will modify the keys of this dictionary. This allows 
# to no have to pass certain usseless keys
# in addition it will ignore keys that they are not in the dictionary, to avoid crashing

k0 = 100.
k1 =  1.
k2 = 1.
k3 = 1.
cellcycleduration = 1.
m0 = 0
p0 = 0

basedictionary["Omega"] = 1
basedictionary["k0"] = k0
basedictionary["k1"] = k1
basedictionary["k2"] = k2
basedictionary["k3"] = k3
basedictionary["T"] = cellcycleduration
basedictionary["m0"] = m0
basedictionary["p0"] = p0
basedictionary["ffwrite"] = 1
basedictionary["timelapse"] = 1
basedictionary["totaltime"] = 100
basedictionary["dt"] = 0.01
basedictionary["runtype"] = 0
basedictionary["SEED"] = -1
basedictionary["runtimes"] = 1
basedictionary["stocycle"] = 0
basedictionary["phasenumber"] = 2
basedictionary["presynthesis"] = 1
basedictionary["cellphaserates"] = 1
basedictionary["mRNAgeo"] = -1

mRNAbeforeDivide = []
mRNAafterDivide = []
protafterDivide = []
protbeforeDivide = []

time = 0
celltime = 0
m = 0
p = 0

######################## Functions to wrap the C code

def CreateInputFile_C(input_dic):

	with open('input_gillespie_cyclo.in', 'w') as inputfile:
		for key in basedictionary:
			if key in input_dic:
				if hasattr(input_dic[key], "__len__"): # in case that is a list, write every element
					inputfile.write(key)
					for el in input_dic[key]:	
						inputfile.write(' '+str(el))
					inputfile.write('\n')
				else:
					inputfile.write(key+' '+str(input_dic[key])+'\n')
			else:
				inputfile.write(key+' '+str(basedictionary[key])+'\n')

def LoadOutputfile_C():

	return np.loadtxt('output_gillespie_cyclo.out')

def RunTime_C(input_dic = basedictionary):

	CreateInputFile_C(input_dic)
	os.system('./cyclo_2states input_gillespie_cyclo.in output_gillespie_cyclo.out > gillespie_cyclo.log')
	return LoadOutputfile_C()

def RunTime_C_Pop(input_dic = basedictionary):

	CreateInputFile_C(input_dic)
	os.system('./cyclo_2states_pop input_gillespie_cyclo.in output_gillespie_cyclo.out > gillespie_cyclo.log')
	return LoadOutputfile_C()


###################################################
