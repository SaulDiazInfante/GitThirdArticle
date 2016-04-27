#===============================================================================
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
#===============================================================================
#
#		Tonatihu data
#
#===============================================================================
import numpy as np
import matplotlib.pylab as plt
U1=np.load("UPathsShortTime.npy")
MeanUShortTime = np.mean(U1, axis=0)
MSUShortTime = np.sqrt(np.mean((U1**2),axis=0))
del U1
V1=np.load("VPathsShortTime.npy")
MeanVShortTime = np.mean(V1, axis=0)
MSVShortTime = np.sqrt(np.mean((V1**2),axis=0))
del V1
np.save("MeanUShortTime.npy",\
		np.transpose(np.array([MeanUShortTime[:], MeanVShortTime[:]]))
)
np.save("MSUShortTime.npy",\
		np.transpose(np.array([MSUShortTime[:], MSVShortTime[:]]))
)
#
plt.figure()
plt.plot(MeanUShortTime)
plt.plot(MSUShortTime)
plt.figure()
plt.plot(MeanVShortTime)
plt.plot(MSVShortTime)
plt.show()
U2=np.load("UPathsLongTimeFiltered.npy")
MeanULongTime = np.mean(U2, axis=0)
MSULongTime = np.sqrt(np.mean((U2**2),axis=0))
del U2

V2=np.load("VPathsLongTimeFiltered.npy")
MeanVLongTime = np.mean(V2, axis=0)
MSVLongTime = np.sqrt(np.mean((V2**2),axis=0))
del V2
np.save("MeanULongTime.npy",\
		np.transpose(np.array([MeanULongTime[:], MeanVLongTime[:]]))
		)

np.save("MSULongTime.npy",\
		np.transpose(np.array([MSULongTime[:], MSVLongTime[:]]))
		)
plt.figure()
plt.plot(MeanULongTime)
plt.plot(MSULongTime)
plt.figure()
plt.plot(MeanVLongTime)
plt.plot(MSVLongTime)
plt.show()
