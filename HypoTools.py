from math import sqrt
from statistics import stdev
from scipy.stats import norm, t, chi2
from numpy import array
    

class sampleGroup():
    def __init__(self, n, xbar, sigma, popStdevKnown):
        '''
sigma --->  Positive real number.
            It stands for sample deviation if popStdevKnown parameter is False,
            stands for population deviation if popStdevKnown parameter is True
popStdevKnown ---> True if population deviation is known, False if population deviation is unknown
'''
        self.n = n
        self.xbar = xbar
        self.sigma = sigma
        self.is_popVarianceKnown = popStdevKnown

def used_test(sampleGroup):
    if sampleGroup.is_popVarianceKnown:
        return("z-score test is applied.")
    else:
        return("t-student test is applied.")

#----------------------------------------------------------------------------------------------------------------------
class sampleProp():
    def __init__(self, prop, n):
        self.prop = prop
        self.n = n
        
#======================================================================================================================
        
def zalpha_RTT(alpha):
    zalpha = norm.ppf(1-alpha)
    return(zalpha)

def zalpha_LTT(alpha):
    zalpha = norm.ppf(alpha)
    return(zalpha)

def zalpha_2TT(alpha):
    zalpha = abs(norm.ppf(alpha/2))
    return(zalpha)

zalphaTails_Dict = {'>': zalpha_RTT,
                    '<': zalpha_LTT,
                    '!=': zalpha_2TT}
#======================================================================================================================

def talpha_RTT(alpha, dof):
    talpha = t.ppf(1-alpha, dof)
    return(talpha)

def talpha_LTT(alpha, dof):
    talpha = t.ppf(alpha, dof)
    return(talpha)

def talpha_2TT(alpha, dof):
    talpha = abs(t.ppf(alpha/2, dof))
    return(talpha)

talphaTails_Dict = {'>': talpha_RTT,
                    '<': talpha_LTT,
                    '!=': talpha_2TT}
#======================================================================================================================

def decn_RTT(vtest, valpha):
    if vtest < valpha: return(False)  #Fail to reject null hypothesis
    if vtest >= valpha: return(True)  #Success to reject null hypothesis

def decn_LTT(vtest, valpha):
    if vtest <= valpha: return(True)    #Success to reject null hypothesis
    if vtest > valpha: return(False)  #Fail to reject null hypothesis

def decn_2TT(vtest, valpha):
    if -valpha<vtest<valpha: return(False)                  #Fail to reject null hypothesis
    if (vtest<-valpha or vtest>valpha): return(True)        #Success to reject null hypothesis

Decision = {'>': decn_RTT,
            '<': decn_LTT,
            '!=': decn_2TT}
#======================================================================================================================
