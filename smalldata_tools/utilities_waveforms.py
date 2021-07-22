import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.fft
import math
from scipy.signal import argrelmax, argrelmin

#sys.path.insert(0, '/cds/home/a/akamalov/MRCO_Scripts/SharedMemoryAnalysis/library')

#calculatePSD takes a time-domain input, dataIn, and calculates the PSD in frequency space.  The frequency components are trimmed to only return non-degenerate positive frequencies.
def calculatePSD(dataIn):
    #calculate the FFT of the raw trace
    discreteFourierTransform = scipy.fft.fft(dataIn)
    #calculate the PSD of the trace by first taking the absolute value, and then squaring the FFT
    absOfFFT = np.absolute(discreteFourierTransform)
    PSD = np.square(absOfFFT)
    #output fourier space components should be trimmed to only include the non-degenerate terms
    PSDOfTrace = cutFourierScaleInHalf(PSD)

    #return the calculated power spectral density
    return PSDOfTrace

#convert the time-axis of the collected trace into an axis of frequencies.
def convertTimeAxisToFrequencies(timeAxis):
    #make sure that the timeAxis is a numpy array
    timeAxisArray = np.asarray(timeAxis)
    #calculate the number of samples and the time differential between each step
    numSamples = len(timeAxisArray)
    differentialSpacing = (timeAxisArray[-1] - timeAxisArray[0])/(numSamples-1)
    #use these values to call fftfreq.  The length of freqAxis is the same as the length of timeAxis.
    freqAxis = np.fft.fftfreq(numSamples, differentialSpacing)

    return freqAxis

#method for cutting the frequency spectrum in half.  Fourier transforms produce identical positive and negative frequency components.  Plotting both is unnecessary, so it is useful to have a method that scales an axis in half.
def cutFourierScaleInHalf(inputFourierSpaceArray):
    #there's a slightly different procedure depending on whether the length of the array is even or odd
    numElements = len(inputFourierSpaceArray)
    if((numElements % 2) == 0):
        #the number of elements is even
        #keep the first half of the frequency-space axis.  the second half is redundant
        halfwayIndex = int(numElements/2)
        #also drop the 0th term (the DC component).  keep all remaining non-degenerate frequency components
        nonDegenerateAxis = inputFourierSpaceArray[1:halfwayIndex]
    else:
        #the number of elements is odd
        #first, find the inflection point of the FFT axis (final index before negative frequencies are included)
        inflectionIndex = int(math.ceil(numElements/2))
        #pick out the section of the fourier space axis to keep.  include all frequencies before the inflection point, but drop the 0th term (the DC term)
        nonDegenerateAxis = inputFourierSpaceArray[1:inflectionIndex]

    return nonDegenerateAxis


#compute and return the convolution of a supplied vector, 'vectorToConvolve', when the vector is convolved with a gaussian filter of length 'filterLength'.  The convolution is done using he mode='same' mode of np.convolve.  The filter is a normal distribution described by a linear spacing that spans from -1*maxSigma to maxSigma, such that the end points are not included in the linear spacing.
def convolveWithGaussian(vectorToConvolve, filterLength, maxSigma=3):
    #create a gaussian of filterLength.  Filter should have a uniform spacing of standard deviation values with boundaries of -minSigma to +maxSigma
    tempStdDevSpacing = np.linspace(-1*maxSigma, maxSigma, filterLength + 2)
    #drop endpoints of linspace
    stdDevSpacing = tempStdDevSpacing[1:]
    #compute a normalized filter
    gaussianFilter = np.exp(-0.5*(np.square(stdDevSpacing)))/(np.sqrt(2*3.14159265359))
    gaussianFilter = gaussianFilter/np.sum(gaussianFilter)
    #calculate and return convoluted vector of 'same' length
    convolvedVector = np.convolve(vectorToConvolve, gaussianFilter, mode='same')
    if vectorToConvolve.size < gaussianFilter.size:
        toCut = gaussianFilter.size - vectorToConvolve.size
        if toCut % 2 == 0:
            convolvedVectorTemp = convolvedVector[int(toCut/2):int(toCut/2 + vectorToConvolve.size)]
        else:
            convolvedVectorTemp = convolvedVector[int(toCut/2 + 0.5):int(toCut/2 + 0.5 + vectorToConvolve.size)]
        convolvedVector = convolvedVectorTemp

    return convolvedVector

#filter out the frequencies of an input below some thresholdFrequency, as provided in Hz.  By default is close to an ideal high pass filter.  The signal to be filtered is provided as timeDependentFunction, sampled with time steps of timeSpacing provided in seconds.  If a non-ideal filter is desired, can provide power=1 for a high pass filter with 20dB/octave loss, power=2 for a 40dB/octave loss, etc.
def filterOutFrequenciesBelowThreshold(timeDependentFunction, timeSpacing, thresholdFrequency, power=25):
    #figure out which frequncies fall outside of the threshold
    signalSize = timeDependentFunction.size
    freqAxis = scipy.fft.fftfreq(signalSize, timeSpacing)
    freqsTransferFunction = complex(0,1)*np.power(freqAxis/thresholdFrequency, power)/(1 + complex(0,1)*np.power(freqAxis/thresholdFrequency, power))
    #create a fourier transform only including the values above threshold
    allFrequencies = scipy.fft.fft(timeDependentFunction)
    onlyValidFreqs = allFrequencies * freqsTransferFunction#multiply all frequencies by array which says whether the frequenices pass the threshold
    #create and return the time-dependent signal based only on frequencies that pass the filter
    filteredTimeDependentFunction = scipy.fft.ifft(onlyValidFreqs)
    #get rid of the imaginary parts.  this is a bit sketchy but I think it's slightly better than returning a physical time dependent value with complex values.
    filteredTimeDependentFunction = np.real(filteredTimeDependentFunction)
    
    return filteredTimeDependentFunction

#filter out the frequencies of an input above some thresholdFrequency, as provided in Hz.  By default is close to an ideal low pass filter.  The signal to be filtered is provided as timeDependentFunction, sampled with time steps of timeSpacing provided in seconds.  If a non-ideal filter is desired, can provide power=1 for a low pass filter with 20dB/octave loss, power=2 for a 40dB/octave loss, etc.
def filterOutFrequenciesAboveThreshold(timeDependentFunction, timeSpacing, thresholdFrequency, power=25):
    #figure out which frequncies fall outside of the threshold
    signalSize = timeDependentFunction.size
    freqAxis = scipy.fft.fftfreq(signalSize, timeSpacing)
    freqsTransferFunction = 1/(1 + complex(0,1)*np.power(freqAxis/thresholdFrequency, power))
    #create a fourier transform only including the values above threshold
    allFrequencies = scipy.fft.fft(timeDependentFunction)
    onlyValidFreqs = allFrequencies * freqsTransferFunction#multiply all frequencies by array which says whether the frequenices pass the threshold
    #create and return the time-dependent signal based only on frequencies that pass the filter
    filteredTimeDependentFunction = scipy.fft.ifft(onlyValidFreqs)
    #get rid of the imaginary parts.  this is a bit sketchy but I think it's slightly better than returning a physical time dependent value with complex values.
    filteredTimeDependentFunction = np.real(filteredTimeDependentFunction)
    
    return filteredTimeDependentFunction

#This method helps filterout artifact frequencies associated with digitizer readouts.  Some digitizers consist of multiple sub-digitizer units that interweave their sampling, but report a different baseline.  The result of this is a high frequency zig-zagging.  This method is designed to help eliminate that zig-zagging by cutting out the spikes seen in the fourier transform.
#signal is the raw trace, in time domain and timeAxis is the associated time vector.  frequenciesToEliminate is the array/list of frequencies (in Hz) that need to be eliminated.  replacementFrequenciesIndexOffsets is the list of indices surrounding the index of the targeted frequency, whose associated frequency values get averaged together to form the replacement value.
def eliminateListedSpikeFrequenciesFromSignal(signal, timeAxis, frequenciesToEliminate, replacementFrequenciesIndexOffsets=np.array([-2, -1, 1, 2])):
    #calculate the FFT of the raw trace
    discreteFourierTransform = scipy.fft.fft(signal)
    #calculate the frequency axis associated with provided time axis
    freqAxis = convertTimeAxisToFrequencies(timeAxis)
    #find the index at which frequencies need to be eliminated, and smooth out the frequencies based on provided index
    for i in frequenciesToEliminate:
        #convert from frequency to eliminate to index that represents the frequency, and index that represents the degenerate inverse
        indexToEliminate = indexOfArrayClosestToValue(freqAxis, i)
        indexInverseToEliminate = indexOfArrayClosestToValue(freqAxis, -1*i)
        #calculate the averaged value that will serve as replacement value for the frequency that is being eliminated
        sumReplacement = 0
        for j in replacementFrequenciesIndexOffsets:
            sumReplacement += discreteFourierTransform[indexToEliminate + j]
        replacement = sumReplacement/replacementFrequenciesIndexOffsets.size
        #replace the i'th frequency to eliminate with the calculated replacement, replace degenerate negative frequency as well
        discreteFourierTransform[indexToEliminate] = replacement
        discreteFourierTransform[indexInverseToEliminate] = np.conj(replacement)
    #calculate the time signal with spike frequencies removed.  returning the value ifft(discreteFourierTransform) itself yields an array with complex arguments.  try returning the real value only and hoping that's close enough
    traceWithoutFrequencies = np.real(scipy.fft.ifft(discreteFourierTransform))#*(np.angle(scipy.fft.ifft(discreteFourierTransform))/np.absolute(np.angle(scipy.fft.ifft(discreteFourierTransform))))
    
    #return value
    return traceWithoutFrequencies


#normalize a 1-D numpy array.  If you're reading this comment and think this method is stupid: fine, but your code is more readable.
def normalizeOneDimArray(arrayIn):
    return arrayIn/np.sum(arrayIn)

#accepts an array, arrayIn, and looks for the index of arrayIn which has a value closest to targetValue.  mostly here to help with code legibility
def indexOfArrayClosestToValue(arrayIn, targetValue):
    return (np.abs(arrayIn - targetValue)).argmin()

#bin a raw histogram into the bin width specified by the user.  return a plot line's y-values to resemble histogram blocks, but be of the same length as the input variable 'rawHistogram'
def calculateBinnedHistogramTrace(rawHistogram, binWidth):
    lenFullTrace = rawHistogram.size
    #calculate the number of complete bins that rawHistogram can be binned into, for given binWidth
    numberCompleteBins = int(np.floor(lenFullTrace/binWidth))
    #reshape as much of the rawHistogram trace as possible
    lengthToBeReshaped = numberCompleteBins*binWidth
    reshapedTraceArray = np.reshape(rawHistogram[0:lengthToBeReshaped], [numberCompleteBins, binWidth])
    #use the reshaped trace to simplify calculation of bins.  sum up along axis 1 to sum across the bin width dimension.  in other words, sum up the components of a single bin with width binWidth.
    sumsOfBins = np.sum(reshapedTraceArray, 1)

    #using the binnedTrace, and the unutilized tail end of rawHistogram, stich together an array of the same dimension as rawHistogram, but with values that represent binned data.
    binnedPortionOfTrace = np.repeat(sumsOfBins, binWidth)#account for the binned portion of the trace.
    #stitch on any part of the trace not used in the binning
    unusedTraceTail = rawHistogram[lengthToBeReshaped:lenFullTrace]
    binnedTrace = np.concatenate((binnedPortionOfTrace, unusedTraceTail))#need argument of method to be a tuple of the two arrays to be stitched together.

    return binnedTrace



#helper method to eliminate the spike frequency artifacts caused by how the digitizer works.
def hsdBaselineFourierEliminate(wf, times):
    #setup to look at fourier spectra of the waveform
    freqAxisFull = convertTimeAxisToFrequencies(times)
    freqMax = np.amax(np.absolute(freqAxisFull))
    fftFull = np.absolute(scipy.fft.fft(wf))
    
    #eliminate listed frequencies (as freqMax/i, where i is listed), from the fourier spectra.  This section eliminates frequencies that are somewhat broader.
    frequencyList = []
    dividerListStrong = [8, 5, 5/2, 5/3, 5/4, 4, 2, -1] #use -1 and not 1.  this is because of how FFT's work for frequency axis with an odd number of elements
    for i in dividerListStrong:
        frequency = freqMax/i
        maxInd = indexOfArrayClosestToValue(freqAxisFull, frequency)
        maxInd = np.argmax(fftFull[maxInd-2:maxInd+2]) + maxInd - 2
        for j in range(-5, 6):
            frequencyList.append(freqAxisFull[maxInd + j])
        frequenciesToEliminate = np.array(frequencyList)
    #call a helper method to eliminate the frequencies listed
    filteredTrace = eliminateListedSpikeFrequenciesFromSignal(wf, times, frequenciesToEliminate, replacementFrequenciesIndexOffsets=np.array([-10, -9, -8, -7, 7, 8, 9, 10]))
    
    #repeat the process for the frequency components that are somewhat weaker/narrower
    frequencyList = []
    dividerListWeak = [64/i for i in range(1, 64)]
    for i in dividerListWeak:
        frequency = freqMax/i
        maxInd = indexOfArrayClosestToValue(freqAxisFull, frequency)
        maxInd = np.argmax(fftFull[maxInd-2:maxInd+2]) + maxInd - 2
        for j in range(-1, 2):
            frequencyList.append(freqAxisFull[maxInd + j])
    frequenciesToEliminate = np.array(frequencyList)
    filteredTrace = eliminateListedSpikeFrequenciesFromSignal(filteredTrace, times, frequenciesToEliminate, replacementFrequenciesIndexOffsets=np.array([-6, -5, -4, -3, 3, 4, 5, 6]))

    #return the filteredTrace
    return filteredTrace






########################################
#copy/pasted methods from commonMethods:
#########################################


# #This method helps filterout artifact frequencies associated with digitizer readouts.  Some digitizers consist of multiple sub-digitizer units that interweave their sampling, but report a different baseline.  The result of this is a high frequency zig-zagging.  This method is designed to help eliminate that zig-zagging by cutting out the spikes seen in the fourier transform.
# #signal is the raw trace, in time domain and timeAxis is the associated time vector.  frequenciesToEliminate is the array/list of frequencies (in Hz) that need to be eliminated.  replacementFrequenciesIndexOffsets is the list of indices surrounding the index of the targeted frequency, whose associated frequency values get averaged together to form the replacement value.
# def eliminateListedSpikeFrequenciesFromSignal(signal, timeAxis, frequenciesToEliminate, replacementFrequenciesIndexOffsets=np.array([-2, -1, 1, 2])):
#     #calculate the FFT of the raw trace
#     discreteFourierTransform = scipy.fft.fft(signal)
#     #calculate the frequency axis associated with provided time axis
#     freqAxis = convertTimeAxisToFrequencies(timeAxis)
#     #find the index at which frequencies need to be eliminated, and smooth out the frequencies based on provided index
#     for i in frequenciesToEliminate:
#         #convert from frequency to eliminate to index that represents the frequency, and index that represents the degenerate inverse
#         indexToEliminate = indexOfArrayClosestToValue(freqAxis, i)
#         indexInverseToEliminate = indexOfArrayClosestToValue(freqAxis, -1*i)
#         #calculate the averaged value that will serve as replacement value for the frequency that is being eliminated
#         sumReplacement = 0
#         for j in replacementFrequenciesIndexOffsets:
#             sumReplacement += discreteFourierTransform[indexToEliminate + j]
#         replacement = sumReplacement/replacementFrequenciesIndexOffsets.size
#         #replace the i'th frequency to eliminate with the calculated replacement, replace degenerate negative frequency as well
#         discreteFourierTransform[indexToEliminate] = replacement
#         discreteFourierTransform[indexInverseToEliminate] = np.conj(replacement)
#     #calculate the time signal with spike frequencies removed.  returning the value ifft(discreteFourierTransform) itself yields an array with complex arguments.  try returning the real value only and hoping that's close enough
#     traceWithoutFrequencies = np.real(scipy.fft.ifft(discreteFourierTransform))#*(np.angle(scipy.fft.ifft(discreteFourierTransform))/np.absolute(np.angle(scipy.fft.ifft(discreteFourierTransform))))
    
#     #return value
#     return traceWithoutFrequencies
    

##################################################
#this is a conventional constant fraction discriminator (CFD) hitfinder.  To tune it, change values for the variables: threshold, CFDOffset, inverseMultiplier.  threshold alters the peak strength above which the algortihm should check for a zero crossing using a conventional CFD method (https://en.wikipedia.org/wiki/Constant_fraction_discriminator).  CFDOffset is the number of waveform bins by which to offset the inverted signal, and inverseMultiplier is the multiplier to set the strength of the inverted signal.
def hitFinder_CFD(dataIn, convFilterLength = 35, CFDOffset = 25, 
                  inverseMultiplier = -0.75, threshold = 4 ):

	#initialize 'hitIndices', which will contain the indices of any hits found in the trace supplied as 'dataIn_Amplitude'
	hitIndices = []

	#normalize the trace to be positive, and have max value of +1
	dataInNormalized = normalizeTrace(dataIn)
	#calculate the variance of the trace
	sigma = np.std(dataInNormalized)

	#calculate an upper threshold above which to look for peaks in the raw trace
	threshold = sigma*threshold
	#return the indices for which the raw data exceeds the threshold.
	dataIn_AboveThreshold_Indices = np.flatnonzero(dataInNormalized > threshold)

	#if it's likely that there are zero hits in this trace, there's no need to perform the remainder of the CFD processing.
	if(len(dataIn_AboveThreshold_Indices) == 0):
		#create an empty array of found hits
		hitIndices = np.asarray([])
		return hitIndices


	#convolve the raw data with a gaussian filter
	convolvedData = convolveWithGaussian(dataInNormalized, convFilterLength)

	#add up an inverse and an offset.  this is the type of approach an electronic CFD performs.
	lengthTrace = len(convolvedData)
	offsetTrace = convolvedData[0:(lengthTrace - CFDOffset)]
	inverseTrace = inverseMultiplier * convolvedData[CFDOffset:lengthTrace]
	#traditional CFD adds a time-offset copy of the trace with an inverser copy of original trace.
	comparedTrace = offsetTrace + inverseTrace
	#shift the region with zero-point crossing to be more centered on the zero cross.  The initial array is found based on being above some amount of standard deviations
	indicesShift = round(CFDOffset * (1 + inverseMultiplier))
	dataIn_AboveThreshold_Indices -= indicesShift

	#call a method which will take the array of indices, and separate that one array into a set of arrays, wherein each array is a continuous set of integers.
	tupleOfRegionIndicesArrays = separateArrayIntoTupleOfContinuousArrays(dataIn_AboveThreshold_Indices)

	#findZeroCrossings for each array of continuous integers
	for ind in range(len(tupleOfRegionIndicesArrays)):
		seriesToProcess = tupleOfRegionIndicesArrays[ind]
		#method 'findZeroCrossings' inspects a series to validate it.  if it's a good zero-crossing, it returns: True, indexOfCrossing.  if it's a bad series, the return is 'False, 0'
		validSeriesFlag, hitIndex = findZeroCrossings(seriesToProcess, comparedTrace)
		#append good hits to the array 'hitIndices'
		if(validSeriesFlag):
			hitIndices.append(hitIndex)

	#there are now a set of found hitIndices.  but these are in respect to the processed comparedTrace.  need to un-shift the indices to represent hits for the actual trace (dataIn_Centered)
	hitIndices = [x + indicesShift for x in hitIndices]

	return hitIndices



#####################################################################################
#support methods for andrei's CFD



#compute and return the convolution of a supplied vector, 'vectorToConvolve', when the vector is convolved with a gaussian filter of length 'filterLength'.  The convolution is done using he mode='same' mode of np.convolve.  The filter is a normal distribution described by a linear spacing that spans from -1*maxSigma to maxSigma, such that the end points are not included in the linear spacing.
def convolveWithGaussian(vectorToConvolve, filterLength, maxSigma=3):
    #create a gaussian of filterLength.  Filter should have a uniform spacing of standard deviation values with boundaries of -minSigma to +maxSigma
    tempStdDevSpacing = np.linspace(-1*maxSigma, maxSigma, filterLength + 2)
    #drop endpoints of linspace
    stdDevSpacing = tempStdDevSpacing[1:]
    #compute a normalized filter
    gaussianFilter = np.exp(-0.5*(np.square(stdDevSpacing)))/(np.sqrt(2*3.14159265359))
    gaussianFilter = gaussianFilter/np.sum(gaussianFilter)
    #calculate and return convoluted vector of 'same' length
    convolvedVector = np.convolve(vectorToConvolve, gaussianFilter, mode='same')
    if vectorToConvolve.size < gaussianFilter.size:
        toCut = gaussianFilter.size - vectorToConvolve.size
        if toCut % 2 == 0:
            convolvedVectorTemp = convolvedVector[int(toCut/2):int(toCut/2 + vectorToConvolve.size)]
        else:
            convolvedVectorTemp = convolvedVector[int(toCut/2 + 0.5):int(toCut/2 + 0.5 + vectorToConvolve.size)]
        convolvedVector = convolvedVectorTemp

    return convolvedVector

#this method is designed to take an array of integers, some of which are continuous, and separate it into a set of arrays wherein each array is a continuous set of integers.  these individual arrays are placed into a tuple that is then returned.
def separateArrayIntoTupleOfContinuousArrays(dataIn_AboveThreshold_Indices):
	#setup the 'first' currentList and the tuple that will be populated
	currentList = []
	tupleOfLists = ()

	#handle the odd case that there is exactly 1 index found.  This is a rarity, but it needs to be handled to avoid error
	if len(dataIn_AboveThreshold_Indices) == 1:
		currentList += dataIn_AboveThreshold_Indices[0]
		tupleOfLists += (currentList,)
	#the cases which matter are the ones that have more than one element, and are handled in the else statement
	else:
		for ind in range(0, len(dataIn_AboveThreshold_Indices) - 1):
			#add the current index to the current list
			currentList.append(dataIn_AboveThreshold_Indices[ind])

			#inspect whether the next element in the list of indices is the start of a new continuous set.  if it is, close out this list
			if (dataIn_AboveThreshold_Indices[ind + 1] - dataIn_AboveThreshold_Indices[ind]) != 1:
				#the next index is a the start of a new continuous set
				tupleOfLists += (currentList,)
				#clear the currentList, so that the next value considered will be the first value in a new array
				currentList = []

		#process the final index in the array, and close out the current list since the list of indices is complete
		currentList.append(dataIn_AboveThreshold_Indices[-1])
		tupleOfLists += (currentList,)

	return tupleOfLists


#method findZeroCrossings inspects the index series in seriesToProcess, and verifies that the associated y-values in comparedTrace are an appropriate rising edge.  if it's a good series, return true and the zero crossing index.  if false, return False and 0
def findZeroCrossings(seriesToProcess, comparedTrace):
	numIndices = len(seriesToProcess)
	if numIndices <= 1:
		#series of length 1 won't have a proper zero crossing and are therefore, not valid zero crossings
		return False, 0
	else:
		#the ideal zero crossing series starts negative and trends positive.  it is good to filter series for validity by verifying this.
		seriesLowest = seriesToProcess[0]
		seriesHighest = seriesToProcess[-1]

		#verify that series crossing isn't too close to either start or end of the trace.  If it is too near to either, can't test whether zero crossing is valid.
		#Note that the condition checks look at seriesLowest - 1 and seriesHighest + 1.  This is because the way the while loops go through below, the loop can cause either indLow or indHigh to go out of bounds of comparedTrace, and then require a call to comparedTrace with an invalid index on the next boolean condition check.
		if (seriesLowest - 1) < 0 or seriesHighest < 0:
			#verify that the seriesToProcess does not include negative integers - that is, that it is not too close to the start of trace to pass the test
			#if it is, return that the series is not valid
			return False, 0
		elif (seriesHighest + 1) >= len(comparedTrace) or seriesLowest >= len(comparedTrace):
			#verify that the seriesToProcess is not too close to the end of the trace - that is, verify it isn't at the cutoff edge of the time axis.
			#if it is, return that the series is not valid
			return False, 0


		#inspect where the series stops being negative
		indLow = seriesLowest
		while (comparedTrace[indLow] < 0) and (indLow <= seriesHighest):
			indLow += 1

		#inspect where the series stops being positive if coming in from the positive side
		indHigh = seriesHighest
		while (comparedTrace[indHigh] > 0) and (indHigh >= seriesLowest):
			indHigh -= 1

		#if indLow and indHigh are adjacent to each other, then the series passed in was a monotonically positive zero-crossing.
		if ((indHigh + 1) == indLow): #the way the while loops are broken out of, it's a valid series if indLow is one value higher than indHigh
			#return true, and the index of the first positive value after the crossing.
			return True, indHigh
		else:
			#this was not a valid series
			return False, 0



#function to normalize a trace to be positive, such that the max of the trace is 1.
def normalizeTrace(dataIn):
	#ensure that the dataIn value is normalized to zero.  Do this by finding a median value of the trace
	dataInNorm = dataIn - np.median(dataIn)
	#normalize the data such that the highest point has absolute value of 1.  First, find the maximal value but also figure out if peak is upwards or downwards going
	maximalValueAbs = np.absolute(np.max(dataInNorm))
	minimalValueAbs = np.absolute(np.min(dataInNorm))
	if(maximalValueAbs > minimalValueAbs):
		#the peak is positive going.  normalize with positive value
		dataInNorm = dataInNorm/maximalValueAbs
	else:
		#the peak is negative going.  normalize with negative value
		dataInNorm = -1*dataInNorm/minimalValueAbs

	return dataInNorm

