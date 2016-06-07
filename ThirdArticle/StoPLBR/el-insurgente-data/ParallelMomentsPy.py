#=======================================================================================================================
# Calculate the mean and mean square of the BR-solution at Short and Long time.
#INPUT: Binary files of paths.
#	U1PathsShortTime.npy
#	U2PathsShortTime.npy
#	U1PathsLongTime.npy
#	U2PathsLongTime.npy
#OUTPUT: Binary Files with the moments paths.
#	MeanUShortTime.npy
#	MSUShortTime.npy
#	MeanULongTime.npy
#	MSULongTime.npy
#=======================================================================================================================

import numpy as np
import matplotlib.pylab as plt
'''
Up1=np.load("U1PathsShortTime.npy")
MeanU1 = np.mean(Up1, axis=0)
MSU1 = np.sqrt(np.mean((Up1**2),axis=0))
del Up1

Up2=np.load("U2PathsShortTime.npy")
MeanU2 = np.mean(Up2, axis=0)
MSU2 = np.sqrt(np.mean((Up2**2),axis=0))
del Up2
'''

Up1=np.load("U1PathsLongTimeFiltered.npy")
MeanU1LongTime = np.mean(Up1, axis=0)
MSU1LongTime = np.sqrt(np.mean((Up1**2),axis=0))
del Up1

Up2=np.load("U2PathsLongTimeFiltered.npy")
MeanU2LongTime = np.mean(Up2, axis=0)
MSU2LongTime = np.sqrt(np.mean((Up2**2),axis=0))
del Up2

np.save("MeanULongTime.npy",\
		np.transpose(np.array([MeanU1LongTime[:], MeanU2LongTime[:]]))
		)

np.save("MSULongTime.npy",\
		np.transpose(np.array([MSU1LongTime[:], MSU2LongTime[:]]))
		)