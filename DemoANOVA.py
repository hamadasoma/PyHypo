from ANOVA import *

group1 = array([8,7,8,6,9,7])
group2 = array([6,7,5,5,8,5])
group3 = array([5,3,3,5,4,3])


samples = [group1, group2, group3]
alpha = 0.1
case = AnovaTest(samples, alpha)
