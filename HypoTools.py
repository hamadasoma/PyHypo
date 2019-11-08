from math import sqrt
from statistics import stdev
from scipy.stats import norm, t, chi2, f
from numpy import array


DecnWords = {True: "SUCCESS to REJECT null hypothesis.",
             False: "FAIL to REJECT null hypothesis."}
    

class SampleGroup():
    def __init__(self, popStdevKnown, n, xbar, sigma):
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

def used_test(SampleGroup):
    if SampleGroup.is_popVarianceKnown:
        return("z-score test is applied.")
    else:
        return("t-student test is applied.")

#----------------------------------------------------------------------------------------------------------------------
class SampleProp():
    def __init__(self, prop, n):
        self.prop = prop
        self.n = n
        
#======================================================================================================================
                               # Normal distributions 
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
                                       # Chi-Square distributions
#======================================================================================================================

def x2alpha_RTT(alpha, dof):
    x2alpha = chi2.ppf(1-alpha, dof)
    return(x2alpha)

def x2alpha_LTT(alpha, dof):
    x2alpha = chi2.ppf(alpha, dof)
    return(x2alpha)

def x2alpha_2TT(alpha, dof):
    x2alpha = alpha/2
    x2alpha_lt = chi2.ppf(x2alpha, dof)
    x2alpha_rt = chi2.ppf(1-x2alpha, dof)
    return((x2alpha_lt, x2alpha_rt))

chi2_alphaTails_Dict = {'>': x2alpha_RTT,
                    '<': x2alpha_LTT,
                    '!=': x2alpha_2TT}

def decn_2TT_chi2(vtest, chi_tuple):
    if(chi_tuple[0] < vtest < chi_tuple[1]): return(False)          #Fail to reject null hypothesis
    if(vtest<chi_tuple[0] or vtest>chi_tuple[1]): return(True)      #Success to reject null hypothesis

Decision_chi2 = {'>': decn_RTT,
                 '<': decn_LTT,
                 '!=': decn_2TT_chi2}
#======================================================================================================================
