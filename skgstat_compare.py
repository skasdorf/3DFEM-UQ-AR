import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
import skgstat as skg
from scipy.optimize import curve_fit
from pprint import pprint
from scipy.stats import norm
from gstools import SRF, Gaussian, krige


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


# sortArgs = np.argsort(xData)
# xData = xData[sortArgs[::1]]
# yDataRCS = yDataRCS[sortArgs[::1]]





# # xData = np.random.normal(10, 1, 201)
# # xData = np.linspace(0, 5, 201)
# # yData = np.random.normal(2, 0.5, 201)
# yData = np.sort(yData)
# xData = np.sort(xData)

# xData = np.linspace(0, 10, 1000)
# yData = norm.pdf(xData, loc = 5, scale = 1)


####___________________________________ gstools check ____________________________________________________________________



# xSamples = np.array([xData[0], xData[500], xData[998]])
# ySamples = np.array([yDataRCS[0], yDataRCS[500], yDataRCS[998]])

# drift_model = Gaussian(dim= 1, var = 0.1, len_scale = 2)
# drift = SRF(drift_model, seed = 101)


# model = Gaussian(dim = 1, var = 0.1, len_scale = 2)
# krig = krige.Universal(model, cond_pos = xSamples, cond_val = ySamples)

# krig(xData)

# ax = krig.plot()
# # ax.scatter(xSamples, ySamples, color='k', zorder=10, label = 'Conditions')
# # ax.plot()






###__________________________________________scikit gstat and my own version below this ________________________________________________________


# xSamples = np.array([xData[0], xData[125], xData[250], xData[500], xData[750], xData[999]])
# ySamples = np.array([yData[0], yData[125], yData[250], yData[500], yData[750], yData[999]])

xSamples = np.array([xData[5], xData[250], xData[500], xData[750], xData[995]])
ySamples = np.array([yDataRCS[5], yDataRCS[250], yDataRCS[500], yDataRCS[750], yDataRCS[995]])
yData = yDataRCS

# xSamples = np.array([xData[0], xData[50], xData[100], xData[150], xData[200], xData[250], xData[300], xData[350], xData[400], xData[450], xData[500], xData[550], xData[600], xData[650], xData[700], xData[750], xData[800], xData[850], xData[900], xData[950], xData[999]])
# ySamples = np.array([yData[0], yData[50], yData[100], yData[150], yData[200], yData[250], yData[300], yData[350], yData[400], yData[450], yData[500], yData[550], yData[600], yData[650], yData[700], yData[750], yData[800], yData[850], yData[900], yData[950], yData[999]])


xSampLen = len(xSamples)
xDataLen = len(xData)

# #_______________________________________________MATLAB Attempt VARIOGRAM______________________________________________________________________

# maxDist = (xSamples[xSampLen-1] - xSamples[0]) / 1.0
# # maxDist = 4.8
# nBins = 10
# binTol = maxDist / nBins

# #calculate distance matrix
# # dMat = np.zeros((int(xSampLen*(xSampLen+1)/2), 4))
# dMat = np.zeros((xSampLen**2,4))

# #index refers to the combination of i-j
# index = 0
# for i in range(xSampLen):
# 	for j in range(xSampLen):
# 		# if (i <= j): #just takes upper triangular part
# 		dMat[index, 0] = i #ySamples[i]
# 		dMat[index, 1] = j #ySamples[j]
# 		dMat[index, 2] = abs(xSamples[i] - xSamples[j])
# 		dMat[index, 3] = (ySamples[i] - ySamples[j])**2
# 		index +=1


# # print(dMat[:,:])	

# edges = np.linspace(0, maxDist, nBins+1)
# edges[len(edges)-1] = maxDist


# #determine counts in each bin as well as the bin index each i-j combination falls into
# # binIdx = np.zeros((int(xSampLen*(xSampLen+1)/2)))
# binIdx = np.zeros((xSampLen**2))
# binCounts = np.zeros((len(edges)-1))
# for i in range(len(edges)-1):
# 	# for j in range(xSampLen**2):  #matrix sweep
# 	for j in range(len(binIdx)):
# 		if (dMat[j, 2] < edges[i+1] and dMat[j, 2] >= edges[i]):
# 			binCounts[i] += 1
# 			binIdx[j] = i
# variogram = np.zeros((len(edges)-1))

# #sum the values for each bin, then divide by the number of values
# for i in range(len(edges)-1):  #bin sweep
# 	# for j in range(xSampLen**2):  #matrix sweep
# 	for j in range(len(binIdx)):
# 		if (binCounts[i] != 0):
# 			if (binIdx[j] == i):
# 				variogram[int(binIdx[j])] += dMat[j,3]
# 				# print("i: " + str(i) + " j: " + str(j) + " binCounts[i]:" + str(binCounts[i]) + " binIdx[i]: " + str(int(binIdx[i])) + " dMat[j,3]:" + str(dMat[j,3]) + " dMat[j,2]: " + str(dMat[j,2]))
# 	# variogram = variogram / (2.0 * binCounts[i])
# for i in range(len(variogram)):
# 	# if (binCounts[i] != 0):
# 	variogram[i] /= (2.0 * binCounts[i])

# bins = edges[1:len(edges)]# + binTol/2.0

# #trim zero values
# newVariogram = []
# newEdges = []
# newBins = []
# for i in range(len(variogram)):
# 	if ~np.isnan(variogram[i]):
# 		newBins.append(bins[i])
# 		newVariogram.append(variogram[i])
# 		newEdges.append(edges[i])
# bins = np.array(newBins)
# # edges = np.array(newEdges)
# variogram = np.array(newVariogram)
# # variogram = v.data()[1]
# ###---------------------------------------------------------------------------------------------------------------------------------------------------------


# # if (variogram[0] > variogram[1]):
# # 	np.delete(variogram)


# #gauss/exp fitting
# # parameters, covariance = curve_fit(expVar, bins, variogram)
# # p = parameters[0]
# # r = parameters[1]
# # n = parameters[2]

# # d = np.linspace(0, maxDist, 100)

# # variogramOld = variogram
# # variogram = expVar(d, p, r, n)
# #------------------

# #lin fitting
# parameters, covariance = curve_fit(linVar, bins, variogram)
# m = parameters[0]
# # b = parameters[1]

# d = np.linspace(0, maxDist, 10)
# variogramOld = variogram
# variogram = linVar(d, m)
####----------------------------------------------------------------------------------------------------------------------------------------------

drift_model = Gaussian(dim=1, var = 0.1, len_scale = 2)
drift = SRF(drift_model, seed=101)

ySamples = ySamples

model = Gaussian(dim=1, var=0.1, len_scale=2)
gsKrig = krige.Universal(model, xSamples, ySamples, "linear")
gsKrig(xData)



v = skg.Variogram(xSamples, ySamples, model = 'spherical')
# v.plot()
ok = skg.OrdinaryKriging(v)
z = ok.transform(xData, np.ones(xData.shape)*1.0)


plt.figure(3)
ax = gsKrig.plot()
ax.scatter(xData, z, marker = 'x', label = 'predicted values')
ax.scatter(xData, yData, marker = '.', label = 'actual values')
plt.legend()


plt.figure(2)
# plt.scatter(bins, variogramOld)
# plt.plot(d, variogram)
data = v.data()
plt.plot(data[0], data[1])
plt.scatter(v.bins, v.experimental)
plt.xlabel("h", fontsize = 18)
plt.ylabel("$\gamma(h)$", fontsize = 18)
plt.show()


#this makes the calculated variogram put into normal method
# variogram = v.data()
# xDistCheck = variogram[0]
# variogram = variogram[1]

#__________________________________ Old Attempt_________________________________________________________
#create variogram
# h = 1e8
# delta = []
# for i in range(xSampLen):
# 	for j in range(xSampLen):
# 		if (abs(xSamples[i] - xSamples[j]) < h and abs(xSamples[i] - xSamples[j]) > 0):
# 			h = abs(xSamples[i] - xSamples[j])

# h = 0.5
# end = h
# hVec = []
# i = 4
# while end < (xSamples[xSampLen-1] - xSamples[0]):
# 	end = h*i

# 	print(end)

# 	hVec.append(end)
# 	i += 1

# # for i in range(xSampLen-1):
# # 	delta.append(abs(h[i+1] - h[i]))

# h = np.array(hVec)
# # h = np.linspace(0, xSamples[xSampLen-1] - xSamples[0], 5)
# variogram = np.zeros((len(h)))
# varioCount = np.zeros((len(h)))

# for i in range(xSampLen):
# 	for j in range(xSampLen):
# 		xDist = abs(xSamples[i] - xSamples[j])
# 		[value, index] = find_nearest(h, xDist)
# 		print("i: " + str(i) + " j: " + str(j) + " index: " + str(index))
# 		variogram[index] += (ySamples[i] - ySamples[j])**2
# 		varioCount[index] += 1

# variogram = np.divide(variogram, (2.0 * varioCount))
# plt.scatter(h, variogram)
# plt.show()



# # # _________________________Kriging System__________________________________________________

# variance = np.zeros((xSampLen+1,xSampLen+1))
# for i in range(xSampLen+1):
# 	for j in range(xSampLen+1):

# 		if i == xSampLen:
# 			variance[i][j] = 1
# 		elif j == xSampLen:
# 			variance[i][j] = 1
# 		else:
# 			#work on calculating variogram as a function of h first.  Then fill in the variance matrix.
# 			# xDist = abs(xSamples[i] - xSamples[j])
# 			xDist = np.sqrt((xSamples[i]-xSamples[j])**2)
# 			# xDist = xDist + 1
# 			[value, index] = find_nearest(d, xDist)
# 			# [value, index] = find_nearest(edges, xDist)
# 			variance[i][j] = variogram[index]


# variance[xSampLen][xSampLen] = 0

# varVec = np.zeros((xSampLen+1, 1))
# iVec = np.zeros((xSampLen+1, 1))
# xDistVec = np.zeros((xSampLen+1, 1))
# xPredicted = []
# for k in range(xDataLen): #each pass in this loops predicts a yVal from xData values
# 	#xData[k] is one of the 1000 monte carlo points

# 	for i in range(xSampLen):  #this creates the cov vec
# 		#xSamples[i] is x_i, xData[k] is x* or one of the 1000 mc points
# 		xDist = np.sqrt((xSamples[i]-xData[k])**2)
# 		# xDist = abs(xSamples[i] - xData[k])
# 		# xDist = np.sqrt(xSamples[i]**2 + xData[k]**2)

# 		# print(xDist)
# 		# [value, index] = find_nearest(edges, xDist)
# 		[value, index] = find_nearest(d, xDist)

# 		iVec[i] = index
# 		varVec[i] = variogram[index]
# 		xDistVec[i]  = xDist
# 		# varVec[i] = 1/2.0*(ySamples[i] - yData[k])**2
# 		# varVec[i] = 1/2.0 * (xSamples[i] - xData[k])**2
# 	varVec[xSampLen] = 1

# 	weights = 0
# 	weights = np.matmul(np.linalg.inv(variance), varVec)
# 	# factor = weights[len(weights)-1]
# 	# if (factor < 0):
# 	# 	factor = -1.0/2.0 * np.sqrt(abs(factor))
# 	# else:
# 	# 	factor = 1.0/2.0 * np.sqrt(abs(factor))
# 	# print(weights)
# 	lagrange = weights[-1]
# 	weights = np.delete(weights, len(weights)-1)

# 	# for i in range(weights.size()):
# 	# 	if weights[i] = 1:

# 	# print(weights)
# 	# print(xSamples)

# 	# xPredicted.append(np.matmul(np.transpose(weights), xSamples))
# 	xPredicted.append(sum(weights*ySamples))
# 	# print("xPredicted: " + str(xPredicted) + "\n")
# 	# print("weights" + str(weights) + "\n")

# plt.figure(4)
# xPredicted = np.array(xPredicted)
# plt.scatter(xData, yData, marker = '.', label = 'actual')
# plt.scatter(xData, xPredicted, label = 'predicted')
# # plt.ylim(0, 10)
# plt.legend()
# plt.show(block=False)