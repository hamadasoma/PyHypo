from statistics import variance
from itertools import combinations
from scipy.stats import f
from numpy import array, append

class AnovaTest():
    '''
Null hypothesis: Population mean for all samples is the same
Alternative hypothesis: At least one sample has a different population mean
'''
    def __init__(self, samples, alpha):
        #samples = [sample1, sample2, ...., samplek]
        #SamplesMeans = array([xbar1, xbar2, ..., xbark])
        #SamplesVariances = array([var1, var2, ..., vark])
        self.k = len(samples)
        self.alpha = alpha
        SamplesMeans = array([])
        SamplesVariances = array([])
        SamplesSizes = array([])
        for asample in samples:
            SamplesMeans = append(SamplesMeans, asample.mean())
            SamplesVariances = append(SamplesVariances, variance(asample, xbar = asample.mean()))
            SamplesSizes = append(SamplesSizes, asample.size)
        self.xbargm = SamplesMeans.mean()
        self.xbars = SamplesMeans
        self.S_squares = SamplesVariances
        self.n = SamplesSizes
        #Falpha calculation
        self.dof1 = self.k - 1
        self.dof2 = self.n.sum() - self.k
        self.Falpha = f.ppf(1-self.alpha, self.dof1, self.dof2)
        #Ftest calculation
        self. SB_square = (self.n * (self.xbars-self.xbargm)**2).sum() / (self.k - 1)
        self. SW_square = ((self.n - 1) * self.S_squares).sum() / (self.n - 1).sum()
        self.Ftest = self. SB_square / self. SW_square
        #Decision --- according to LTT of F-distribution curve.
        if self.Ftest >= self.Falpha : self.reject_H0 = True    #SUCCESS to reject null hypothesis
        else : self.reject_H0 = False                           #FAIL to reject null hypothesis

    def ScheffeTest(self):        
        if self.reject_H0 == False: return("Population mean for all samples is the same. No need for test")
        else:
            #Calculating Fschalpha
            self.Fschalpha = self.Falpha * (self.k - 1)
            #Calculating Fschtest
            comaprison = {'Fschalpha': self.Fschalpha} 
            for pair in combinations(range(1, self.k + 1), 2):
                numer = (self.xbars[pair[0]-1] - self.xbars[pair[1]-1])**2
                denom = (1/self.n[pair[0]-1]) + (1/self.n[pair[1]-1])
                Fschtest = (numer / denom) / self.SW_square
                comaprison[pair] = (Fschtest, Fschtest<self.Fschalpha)
                #True---> Within the same population mean.    #False---> Different population mean. 
            return(comaprison)
            
        
        
            

