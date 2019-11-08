from HypoTools import *

'''
   Hypothesis test for means for two populations:
   1 - Population with known variance.
   2 - Population with unknown variance.
   '''
'''
   Two populations match
   1 - Population with unknown variance.
   '''

'''
   Hypothesis test for proportions for two populations
   '''

'''
   Hypothesis test for variances difference for two populations
   '''
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

'''
   Hypothesis test for means for two populations:
   1 - Population with known variance.
   2 - Population with unknown variance.
   '''

class TwoPopulations_MeansEstimation():
    def __init__(self, SampleGroup1, SampleGroup2, alpha, H1_oper, D=0):
        '''
SampleGroup1, SampleGroup2 ---> Two variable cutomized object class of type SampleGroup. Refer to HypoTools.py for details.
alfa ---> Level of significance.
H1_oper ---> Alternative Hypothesis operator: '>', '<', or '!='
D ---> meu1 - meu2, normally zero
'''
        self.SampleGroup1 = SampleGroup1
        self.SampleGroup2 = SampleGroup2
        self.alpha = alpha
        self.H1_oper = H1_oper
        self.D = D
        self.used_test = used_test(SampleGroup1)
        xbar1 = self.SampleGroup1.xbar
        xbar2 = self.SampleGroup2.xbar
        n1 = self.SampleGroup1.n
        n2 = self.SampleGroup2.n
        if self.SampleGroup1.is_popVarianceKnown:     #Population variance is known
            #Calculate zalpha
            self.zalpha = zalphaTails_Dict[self.H1_oper](self.alpha)
            #Calculate ztest            
            sigma1 = self.SampleGroup1.sigma
            sigma2 = self.SampleGroup2.sigma            
            self.ztest = ((xbar1-xbar2) - self.D)/sqrt(sigma1**2/n1 + sigma2**2/n2)
            #Decision
            decn = Decision[self.H1_oper](self.ztest, self.zalpha)
            if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
            else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.          
        else:                                        #Population variance is not known
            S1 = self.SampleGroup1.sigma
            S2 = self.SampleGroup2.sigma
            #Calculate degree of freedom
            dof_numer = (S1**2/n1 + S2**2/n2)**2
            dof_denom = (1/(n1-1))*(S1**2/n1)**2 + (1/(n2-1))*(S2**2/n2)**2
            self.dof = dof_numer/dof_denom            
            #Calculate talpha
            self.talpha = talphaTails_Dict[self.H1_oper](self.alpha, self.dof)            
            #Calculate ttest
            self.ttest = ((xbar1-xbar2) - self.D)/sqrt(S1**2/n1 + S2**2/n2)
            #Decision
            decn = Decision[self.H1_oper](self.ttest, self.talpha)
            if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
            else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.        
        
    def reporter(self):
       report = {'Population variance known': self.SampleGroup1.is_popVarianceKnown,
                 'alpha': self.alpha}
       if self.SampleGroup1.is_popVarianceKnown:     #Population variance is known
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
        

#======================================================================================================

'''
   Two populations match
   1 - Population with unknown variance.
       
+------+-----------+-----------+--------+
  i    | outcome1  |  outcome2 | di
+------+-----------+-----------+--------+
  1    |  x11      |  x21      |x11-x21
+------+-----------+-----------+--------+
  2    |  x12      |  x22      |x12-x22
+------+-----------+-----------+--------+
  3    |  x13      |  x23      |x13-x23
+------+-----------+-----------+--------+
+------+-----------+-----------+--------+
+------+-----------+-----------+--------+
  n    |  x1n      |  x2n      |x1n-x2n
+------+-----------+-----------+--------+

'''  

class TwoPopulations_OutcomesDiffsEstimation():

    def __init__(self, outcome1, outcome2, alpha, H1_oper, D=0, sigmad='na'):
        '''
outcome1, outcome2 ---> type: numpy arrays of equal sizes.
sigmad:  'na';(default value)---> unknown population deviation, number--->known population deviation
'''
        self.outcome1 = outcome1
        self.outcome2 = outcome2
        self.n = outcome1.size        
        self.H1_oper = H1_oper
        self.alpha = alpha
        self.D = D
        self.sigmad = sigmad        
        if self.sigmad != 'na':   #POPULATION DEVIATION KNOWN
            Sd = self.sigmad / sqrt(self.n)
            pass #z-score
        if self.sigmad == 'na':  #POPULATION DEVIATION NOT KNOWN
            #Calculate ttest
            di = self.outcome1 - self.outcome2
            dbar = di.mean()
            Sd = stdev(di, dbar)
            self.ttest = (dbar-self.D) / (Sd/sqrt(self.n))            
            #Calculate talpha
            self.dof = self.n - 1            
            self.talpha = talphaTails_Dict[self.H1_oper](self.alpha, self.dof)            
            #Decision
            decn = Decision[self.H1_oper](self.ttest, self.talpha)
            if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
            else: self.reject_H0 = False     #FAIL to REJECT null hypothesis. 
    
    def reporter(self):
        report = {'alpha': self.alpha,
                  'ttest': self.ttest,
                  'DoF': self.dof,
                  'talpha': self.talpha,
                  'Test result': DecnWords[self.reject_H0]}       
        
        for i in report.keys():
            print(f"{i}  :  ", report[i])
            
#===============================================================================================
'''
   Hypothesis test for proportions for two populations
   '''

class TwoPopulations_ProportionDiffsEstimation():
    def __init__(self, SampleProp1, SampleProp2, alpha, H1_oper):
        '''
SampleProp1, SampleProp2 ---> Two variable cutomized object class of type SampleProp. Refer to HypoTools.py for details.
'''
        self.prop1 = SampleProp1.prop
        self.n1 = SampleProp1.n
        self.prop2 = SampleProp2.prop
        self.n2 = SampleProp2.n
        self.alpha = alpha
        self.H1_oper = H1_oper
        #Calculate zalpha
        self.zalpha = zalphaTails_Dict[self.H1_oper](self.alpha)        
        #Calculate ztest
        prop1 = self.prop1
        n1 = self.n1
        prop2 = self.prop2
        n2 = self.n2
        poolProp = (prop1*n1 + prop2*n2)/(n1+n2)
        ztest_numer = prop1 - prop2
        ztest_denom = sqrt(poolProp*(1-poolProp)*(1/n1 + 1/n2))
        self.poolProp = poolProp
        self.ztest = ztest_numer / ztest_denom        
        #Decision
        decn = Decision[self.H1_oper](self.ztest, self.zalpha)
        if decn: self.reject_H0 = True   #SUCCESS to REJECT null hypothesis.
        else: self.reject_H0 = False     #FAIL to REJECT null hypothesis.   

    def reporter(self):
        report = {'alpha': self.alpha,
                  'Pool proportion': self.poolProp, 
                  'zalpha': self.zalpha,
                  'ztest': self.ztest,
                  'Test result': DecnWords[self.reject_H0]}
        for i in report.keys():
            print(f"{i}  :  ", report[i])


#=============================================================================================
'''
   Hypothesis test for variances difference for two populations
   '''
class TwoPopulations_VariancesDiffEstimation():
    '''
H0:  Pop variance 1 = Pop variance 2
H1:  Pop variance 1 != Pop variance 2
'''
    def __init__(self, xSampleGroup1, xSampleGroup2, alpha):
        self.alpha = alpha
        #Calcualting ftest
        Sa = max(xSampleGroup1.sigma, xSampleGroup2.sigma)
        Sb = min(xSampleGroup1.sigma, xSampleGroup2.sigma)        
        self.ftest = Sa**2/Sb**2
        #Calcualting falpha
        if xSampleGroup1.sigma > xSampleGroup2.sigma :
            dofa = xSampleGroup1.n - 1
            dofb = xSampleGroup2.n - 1
        else:
            dofa = xSampleGroup2.n - 1
            dofb = xSampleGroup1.n - 1
        self.falpha = f.ppf(1-alpha, dofa, dofb)
        #Dicision   --- Always RTT
        if self.ftest>self.falpha:
            self.reject_H0 = True   #SUCCESS to REJECT null hypothesis ---> variances are not equal
        else:
            self.reject_H0 = False  #FAIL to REJECT null hypothesis   ---> variances are equal

    def reporter(self):
        report = {'alpha': self.alpha,
                  'ftest': self.ftest, 
                  'falpha': self.falpha,                  
                  'Test result': DecnWords[self.reject_H0]}
        for i in report.keys():
            print(f"{i}  :  ", report[i])
            
        

