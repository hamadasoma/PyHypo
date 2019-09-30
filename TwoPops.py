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
#----------------------------------------------------------------------------------------------------------------------

'''
   Hypothesis test for means for two populations:
   1 - Population with known variance.
   2 - Population with unknown variance.
   '''

class TwoPopulations_MeansEstimation():
    def __init__(self, sampleGroup1, sampleGroup2, alpha, H1_oper, D=0):
        '''
alfa ---> Level of significance.
H1_oper ---> Alternative Hypothesis operator: '>', '<', or '!='
D ---> meu1 - meu2, normally zero
'''
        self.sampleGroup1 = sampleGroup1
        self.sampleGroup2 = sampleGroup2
        self.alpha = alpha
        self.H1_oper = H1_oper
        self.D = D
        self.used_test = used_test(sampleGroup1)
        self.report = {'Report state': "Test not applied yet."}

    def reject_H0(self):
        if self.sampleGroup1.is_popVarianceKnown:     #Population variance is known
            #Calculate zalpha
            zalpha = zalphaTails_Dict[self.H1_oper](self.alpha)
            self.report['zalpha'] = zalpha
            #Calculate ztest
            xbar1 = self.sampleGroup1.xbar
            xbar2 = self.sampleGroup2.xbar            
            sigma1 = self.sampleGroup1.sigma
            sigma2 = self.sampleGroup2.sigma
            n1 = self.sampleGroup1.n
            n2 = self.sampleGroup2.n
            ztest = ((xbar1-xbar2) - self.D)/sqrt(sigma1**2/n1 + sigma2**2/n2)
            self.report['ztest'] = ztest
            #Decision
            decn = Decision[self.H1_oper](ztest, zalpha)
            if decn: self.report['Test result'] = "SUCCESS to REJECT null hypothesis."
            else: self.report['Test result'] = "FAIL to REJECT null hypothesis."
            return decn
        else:                                        #Population variance is not known
            xbar1 = self.sampleGroup1.xbar
            xbar2 = self.sampleGroup2.xbar            
            S1 = self.sampleGroup1.sigma
            S2 = self.sampleGroup2.sigma
            n1 = self.sampleGroup1.n
            n2 = self.sampleGroup2.n
            #Calculate degree of freedom
            dof_domin = (S1**2/n1 + S2**2/n2)**2
            dof_nomin = (1/(n1-1))*(S1**2/n1)**2 + (1/(n2-1))*(S2**2/n2)**2
            dof = dof_domin/dof_nomin
            self.report['Degree of Freedom'] = dof
            #Calculate talpha
            talpha = talphaTails_Dict[self.H1_oper](self.alpha, dof)
            self.report['talpha'] = talpha
            #Calculate ttest
            ttest = ((xbar1-xbar2) - self.D)/sqrt(S1**2/n1 + S2**2/n2)
            self.report['ttest'] = ttest
            #Decision
            decn = Decision[self.H1_oper](ttest, talpha)
            if decn: self.report['Test result'] = "SUCCESS to REJECT null hypothesis."
            else: self.report['Test result'] = "FAIL to REJECT null hypothesis."
            return decn

    def reporter(self):
        self.report ['Report state'] = "Test completed."
        self.report ['Population variance known'] = self.sampleGroup1.is_popVarianceKnown
        self.report['alpha'] = self.alpha
        for i in self.report.keys():
            print(f"{i} :  ", self.report[i])

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
        #sigmad:  'na'---> unknown population deviation, number--->known population deviation
        self.outcome1 = outcome1
        self.outcome2 = outcome2
        self.n = outcome1.size        
        self.H1_oper = H1_oper
        self.alpha = alpha
        self.D = D
        self.sigmad = sigmad
        self.report = {'Report state': "Test not applied yet."}

    def reject_H0(self):
        if self.sigmad != 'na':   #POPULATION DEVIATION KNOWN
            Sd = self.sigmad / sqrt(self.n)
            pass #z-score
        if self.sigmad == 'na':  #POPULATION DEVIATION NOT KNOWN
            #Calculate ttest
            di = self.outcome1 - self.outcome2
            dbar = di.mean()
            Sd = stdev(di, dbar)
            ttest = (dbar-self.D) / (Sd/sqrt(self.n))
            self.report['ttest'] = ttest
            #Calculate talpha
            dof = self.n - 1            
            talpha = talphaTails_Dict[self.H1_oper](self.alpha, dof)
            self.report['talpha'] = talpha
            #Decision
            decn = Decision[self.H1_oper](ttest, talpha)
            if decn: self.report['Test result'] = "SUCCESS to REJECT null hypothesis."
            else: self.report['Test result'] = "FAIL to REJECT null hypothesis."
            return decn

    def reporter(self):
        self.report ['Report state'] = "Test completed."        
        self.report['alpha'] = self.alpha
        for i in self.report.keys():
            print(f"{i} :  ", self.report[i])
            
#===============================================================================================
'''
   Hypothesis test for proportions for two populations
   '''

class TwoPopulations_ProportionDiffsEstimation():
    def __init__(self, sampleProp1, sampleProp2, alpha, H1_oper):
        self.prop1 = sampleProp1.prop
        self.n1 = sampleProp1.n
        self.prop2 = sampleProp2.prop
        self.n2 = sampleProp2.n
        self.alpha = alpha
        self.H1_oper = H1_oper
        self.report = {'Report state': "Test not applied yet."}

    def reject_H0(self):
        #Calculate zalpha
        zalpha = zalphaTails_Dict[self.H1_oper](self.alpha)
        self.report['zalpha'] = zalpha
        #Calculate ztest
        prop1 = self.prop1
        n1 = self.n1
        prop2 = self.prop2
        n2 = self.n2
        poolProp = (prop1*n1 + prop2*n2)/(n1+n2)
        ztest_domin = prop1 - prop2
        ztest_nomin = sqrt(poolProp*(1-poolProp)*(1/n1 + 1/n2))
        ztest = ztest_domin / ztest_nomin
        self.report['ztest'] = ztest
        #Decision
        decn = Decision[self.H1_oper](ztest, zalpha)
        if decn: self.report['Test result'] = "SUCCESS to REJECT null hypothesis."
        else: self.report['Test result'] = "FAIL to REJECT null hypothesis."
        return decn

    def reporter(self):
        self.report ['Report state'] = "Test completed."        
        self.report['alpha'] = self.alpha
        for i in self.report.keys():
            print(f"{i} :  ", self.report[i])

     
        
        
            
        
        
  
        
        


        
        
        
            
            
            
            
        
