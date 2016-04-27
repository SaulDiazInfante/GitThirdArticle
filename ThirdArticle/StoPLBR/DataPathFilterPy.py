#===============================================================================
#This code filter the data respect to (OC) and (OB) levels.
#INPUT:
#	U1PahtsLongTime.npy,
#	U2PahtsLongTime.npy,
#	binary files with np-arrays of paths.
# OUTPUT:
#		U1PathsLongTimeFiltered.npy,
#		U2PathsLongTimeFiltered.npy,
#	binary files with the filtered data
##===============================================================================

import numpy as np
import matplotlib.pyplot as plt

Up1 = np.load("U1PathsLongTime.npy")
K = Up1.shape[0]
Index= np.arange(K)
TrashU1=[]
for i in Index:
	if Up1[i,:].max()>= 25:
		TrashU1.append(i)
TrashU1Arr = np.array(TrashU1)	
print 'U1 trash paths: '+str(TrashU1Arr.shape[0]) + ' of '+ str(Up1.shape[0])
del Up1
Up2 = np.load("U2PathsLongTime.npy")
TrashU2=[]
for i in Index:
	if Up2[i,:].max()>= 5000:
		TrashU2.append(i)
TrashU2Arr=np.array(TrashU2)
print 'U2 trash paths: '+str(TrashU2Arr.shape[0]) + ' of '+ str(Up2.shape[0])
TrashIndex = np.union1d(TrashU1Arr,TrashU2Arr)
FilterIndex = np.setdiff1d(Index, TrashIndex)
FilterDataUp2=Up2[FilterIndex]
np.save("U2PathsLongTimeFiltered.npy", FilterDataUp2)
print 'Acceptable paths: '+str(FilterDataUp2.shape[0])
print 'Proportion of rejections paths: '+ str(100.0 -100 * np.float(FilterDataUp2.shape[0])/np.float(Up2.shape[0]))
del Up2
del FilterDataUp2
Up1 = np.load("U1PathsLongTime.npy")
FilterDataUp1=Up1[FilterIndex]
np.save("U1PathsLongTimeFiltered.npy", FilterDataUp1)
del Up1
del FilterDataUp1
