import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
import skgstat as skg
from scipy.optimize import curve_fit
from pprint import pprint
from scipy.stats import norm
from gstools import SRF, Gaussian, krige
import gstools as gs


def printArray(arr):
	for i in range(arr.shape[0]): #rows
		for j in range(arr.shape[1]): #columns
			print('%.4f'%arr[i,j], end = ',')
		print('\n')
		

def gaussVar(x, p, r, n):
	y = p*(1 - np.exp(-x**2/(4.0/7.0*r)**2)) + n
	return y


def expVar(x, p, r, n):
	y = p*(1-np.exp(-x/(r/3))) + n
	return y

def sphereVar(x, p, r, n):
	if (x.any() <= r):
		y = p*(3*x/(2*r) - x**3/(3*r**3))+n
	else:
		y = p+n
	return y

def linVar(x, m):
	y = m*x
	return y


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx


def find_min(array, value):
	array = np.asarray(array)

	for i in range(len(array)):
		if ((value - array[i]) >= 0):
			arrVal = array[i]
			iVal = i
	return arrVal, iVal

def readInput(file):
    db1 = []
    db2 = []
    xData = []
    with open(file, 'r') as f:
        for line in f:
            points = list(map(float, line.split()))
            # print(points[0])
            xData.append(float(points[0]))
            db1.append(points[2])
            db2.append(points[3])
    db1 = np.array(db1)
    db2 = np.array(db2)
    mag = np.sqrt(db1**2 + db2**2)/(4.0*np.pi)
    xData = np.array(xData)
    # print(db)


    return xData, mag

def readMonteCarlo(file):
	inputReal = []
	inputImag = []
	valueReal = []
	valueImag = []
	gradReal = []
	gradImag = []
	with open(file, 'r') as f:
		for line in f:
			points = list(map(float, line.split()))
	        # print(points[0])
			inputReal.append(float(points[0]))
			inputImag.append(float(points[1]))
			valueReal.append(float(points[2]))
			valueImag.append(float(points[3]))
			gradReal.append(float(points[4]))
			gradImag.append(float(points[5]))
	inputReal = np.array(inputReal)
	inputImag = np.array(inputImag)
	valueReal = np.array(valueReal)
	valueImag = np.array(valueImag)
	gradReal = np.array(gradReal)
	gradImag = np.array(gradImag)
	return inputReal, valueReal, valueImag


# xData, yData = readInput('../3DFEM-UQ/cpp/ioFiles/output/qoi_monte_carlo_1000_mat.txt')
# yDataRCS = yData

xData, yDataReal, yDataImag = readMonteCarlo("./cpp/ioFiles/output/MonteCarlo/MC3.txt")
yDataRCS = np.sqrt((yDataReal**2 + yDataImag**2) / (4.0 * np.pi))


sortArgs = np.argsort(xData)
xData = xData[sortArgs[::1]]
yDataRCS = yDataRCS[sortArgs[::1]]
yData = yDataRCS

###_________________________________________________________________________________________________


# xSamples = np.array([xData[0], xData[125], xData[250], xData[500], xData[750], xData[999]])
# ySamples = np.array([yData[0], yData[125], yData[250], yData[500], yData[750], yData[999]])

# xSamples = np.array([xData[0], xData[250], xData[500], xData[750], xData[998]])
# ySamples = np.array([yDataRCS[0], yData[250], yDataRCS[500], yData[750], yDataRCS[998]])
# xSamples = np.array([3.5, 3.7, 4.0, 4.2, 4.5, 4.7, 5.0, 5.2, 5.5, 5.7, 6.0, 6.2, 6.5, 6.7, 7.0, 7.2, 7.5, 7.7, 8.0, 8.2, 8.5, 8.7, 9.0, 9.2, 9.5, 9.7, 10])

xSamples = [3.5]
length = 900
for i in range(length):
	xSamples.append((i+1)*6.5 / length + 3.5)
xSamples = np.array(xSamples)

# xSamples = np.array([3.5, 6, 8, 10])
xSampLen = len(xSamples)
ySamples = np.zeros(xSampLen)

for i in range(xSampLen):
	idx = np.abs([xData - xSamples[i]]).argmin()
	ySamples[i] = yDataRCS[idx]


# xSamples = np.array([xData[0], xData[50], xData[100], xData[150], xData[200], xData[250], xData[300], xData[350], xData[400], xData[450], xData[500], xData[550], xData[600], xData[650], xData[700], xData[750], xData[800], xData[850], xData[900], xData[950], xData[999]])
# ySamples = np.array([yData[0], yData[50], yData[100], yData[150], yData[200], yData[250], yData[300], yData[350], yData[400], yData[450], yData[500], yData[550], yData[600], yData[650], yData[700], yData[750], yData[800], yData[850], yData[900], yData[950], yData[999]])



xDataLen = len(xData)


####----------------------------------------------------------------------------------------------------------------------------------------------
#Universal Kriging
# drift_model = Gaussian(dim=1, var = 0.1, len_scale = 2)
# drift = SRF(drift_model, seed=101)

lenScale = (xSamples[xSampLen-1] - xSamples[0])


model = Gaussian(dim=1, var=0.01, len_scale=lenScale)
gsKrig = krige.Universal(model, xSamples, ySamples, "quadratic")
gsKrig(xData)

ax = gsKrig.plot()
# ax.scatter(xData, z, marker = 'x', label = 'predicted values')
ax.scatter(xData, yData, marker = '.', label = 'actual values')
plt.xlabel("Re{$\epsilon_r$}")
plt.ylabel("Monostatic RCS")
plt.legend()


#variogram
# bins = np.array(40)
# model = Gaussian(dim=2, var=2, len_scale=8)
# srf = SRF(model, mean=0, seed=19970221)
# field = srf((xSamples, ySamples))
# bin_center, gamma = gs.vario_estimate((xSamples, ySamples), field, bins)

# ax1 = fit_model.plot(x_max = 40)
# ax1.scatter(bin_center, gamma)


#this makes the calculated variogram put into normal method
# variogram = v.data()
# xDistCheck = variogram[0]
# variogram = variogram[1]
