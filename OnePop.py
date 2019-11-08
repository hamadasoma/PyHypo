from HypoTools import *

'''
One poulation estimation for:
1 - population mean
    i. known population variance
    ii. unknown population variance
--------------------------------------
2 - population proportion
--------------------------------------
3 - population variance
'''
#=======================================================================


'''
One poulation estimation for:
1 - population mean
    i. known population variance
    ii. unknown population variance
    '''
class OnePopulation_MeanEstimation():
    def __init__(self, xSampleGroup, alpha, H1_oper, meu):        
        self.meu = meu
        self.xbar = xSampleGroup.xbar
        self.sigma = xSampleGroup.sigma
        self.n = xSampleGroup.n
        self.is_popVarianceKnown = xSampleGroup.is_popVarianceKnown
        self.alpha = alpha
        self.H1_oper = H1_oper
        self.used_test = used_test(xSampleGroup)
        if (self.is_popVarianceKnown and self.n>30): #Population variance is known
            #Calculate zalpha
            self.zalpha = zalphaTails_Dict[self.H1_oper](self.alpha)
            #Calculate ztest
            self.ztest = (self.xbar - self.meu)/(self.sigma / sqrt(self.n))
            #Decision
            decn = Decision[self.H1_oper](self.ztest, self.zalpha)
            if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
            else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.
        else:                                       #Population variance is not known
            self.dof = self.n - 1
            #Calculate talpha
            self.talpha = talphaTails_Dict[self.H1_oper](self.alpha, self.dof)            
            #Calculate ttest
            self.ttest = (self.xbar - self.meu)/(self.sigma / sqrt(self.n))
            #Decision
            decn = Decision[self.H1_oper](self.ttest, self.talpha)
            if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
            else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.

    def reporter(self):
       report = {'Population variance known': self.is_popVarianceKnown,
                 'alpha': self.alpha}
       if self.is_popVarianceKnown:     #Population variance is known
            report['zalpha'] = self.zalpha
            report['ztest'] = self.ztest
            report['Test result'] = DecnWords[self.reject_H0]
       else:                                         #Population variance is not known
            report['DoF'] = self.dof
            report['talpha'] = self.talpha
            report['ttest'] = self.ttest
            report['Test result'] = DecnWords[self.reject_H0]
       for i in report.keys():
            print(f"{i}  :  ", report[i])

#=======================================================================
'''
--------------------------------------
2 - population proportion
--------------------------------------
'''
class OnePopulation_ProportionEstimation():
    def __init__(self, xSampleProp, alpha, poprop, H1_oper):
        self.prop = xSampleProp.prop
        self.n = xSampleProp.n
        self.poprop = poprop
        self.alpha = alpha
        self.H1_oper = H1_oper
        #Calculate zalpha
        self.zalpha = zalphaTails_Dict[self.H1_oper](self.alpha)
        #Calculate ztest
        self.ztest = (self.prop-self.poprop)/sqrt(self.poprop*(1-self.poprop)/self.n)
        #Decision
        decn = Decision[self.H1_oper](self.ztest, self.zalpha)
        if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
        else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.

    def reporter(self):
        report = {'alpha': self.alpha,
                  'Sample proportion': self.prop,
                  'Tested population proportion': self.poprop,
                  'zalpha': self.zalpha,
                  'ztest': self.ztest,
                  'Test result': DecnWords[self.reject_H0]}
        for i in report.keys():
            print(f"{i}  :  ", report[i])

#==========================================================================

'''
----------------------------------------
3 - population variance
----------------------------------------
'''
class OnePopulation_VarianceEstimation():
    def __init__(self, xSampleGroup, alpha, H1_oper, popSigma):
        self.S = xSampleGroup.sigma
        self.n = xSampleGroup.n
        self.dof = self.n - 1
        self.alpha = alpha
        self.H1_oper = H1_oper
        self.popSigma = popSigma
        #Calculate X2test
        self.X2test = (self.dof)*(self.S/self.popSigma)**2
        #Calculate X2alpha
        self.X2alpha = chi2_alphaTails_Dict[self.H1_oper](self.alpha, self.dof)
        #Decision
        decn = Decision_chi2[self.H1_oper](self.X2test, self.X2alpha)
        if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
        else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.

    def reporter(self):
        report = {'alpha': self.alpha,
                  'X2alpha': self.X2alpha,
                  'X2test': self.X2test,
                  'Test result': DecnWords[self.reject_H0]}
        for i in report.keys():
            print(f"{i}  :  ", report[i])
        

            
