#simple one-way ANOVA analysis
#pass a dictionary of groups, this should allow a nice summary with the means labeled.
#Also want to be able to pass a coresponding dictionary of contrasts

#Test data

import numpy as np
from scipy import stats

groups = {'a1':[5.0, 9.0, 7.0],
 		'a2':  [12.0,10.0,10.0,8.0,11.0,9.0],
 		'a3':  [3.0,6.0,5.0,6.0]}
 		
 
groups = {'a1':[5.0 , 9.0,  7.0,  4.3, 3.2],
 		'a2':  [12.0, 10.0, 10.0, 8.0, 11.0],
 		'a3':  [3.0,  6.0,  5.0,  6.0, 3.2]}

contrasts = {'a2_vs_a1a3':{'a1':1.0,'a2':-2.0,'a3':1.0},'a1_vs_a3':{'a1':-1.0,'a2':0.0,'a3':1.0}}

class ANOVA:
	def __init__(self,groups = groups,contrasts = contrasts):
		groupnames = groups.keys()
		self.sums = dict()
		self.means = dict()
		self.raw_ss = dict()
		self.sss = dict()
		self.dfs = dict()
		self.sds = dict()
		self.sems = dict()
		self.grandmean = 0
		self.grandn = 0
		self.contrasts = contrasts
		
		for name in groupnames:
			self.grandmean += np.sum(groups[name])
			self.grandn += len(groups[name])
			
		self.grandmean /= self.grandn
		self.grandss = 0
		for name in groupnames:
			n = len(groups[name])
			sum = np.sum(groups[name])
			mean = np.sum(groups[name])/n
			sumsq = np.sum([(x-mean)**2 for x in groups[name]])
			sd = np.sqrt(sumsq/(n-1))
			sem = sd/np.sqrt(n)
			
			self.sums.update({name:np.sum(groups[name])})
			self.means.update({name:mean})
			self.sss.update({name:sumsq})
			self.dfs.update({name:n-1})
			self.sds.update({name:sd})
			self.sems.update({name:sem})
			self.grandss += np.sum([(x-self.grandmean)**2 for x in groups[name]])
			
		self.groupss = np.sum([len(groups[name])*(self.means[name]-self.grandmean)**2 for name in groupnames])
		self.subjectss = np.sum([self.sss[name] for name in groupnames])
		self.groupdf = len(groupnames)-1
		self.subjectdf = np.sum([self.dfs[name] for name in groupnames])
		self.subjectms = self.subjectss/self.subjectdf
		self.groupms = self.groupss/self.groupdf
		self.F = self.groupms/self.subjectms
		self.p = stats.f.sf(self.F,len(groups)-1,self.subjectdf)
		self.contrasts_phi = dict()
		self.contrasts_ms = dict()
		self.contrasts_F = dict()
		self.contrasts_p = dict()
		
		for contrast_name in contrasts.keys():
			contrast = contrasts[contrast_name]
			phi = np.sum([self.means[group]*contrast[group] for group in groupnames])**2
			denom = np.sum([(contrast[group]**2)/len(groups[group]) for group in groupnames])
			contrast_ms = phi/denom
			self.contrasts_phi.update({contrast_name:phi})
			self.contrasts_ms.update({contrast_name:contrast_ms})
			self.contrasts_F.update({contrast_name:contrast_ms/self.subjectms})
			self.contrasts_p.update({contrast_name:stats.f.sf(contrast_ms/self.subjectms,1,self.subjectdf)})
	
	def __str__(self):
		retstr = ""
		retstr += "Omnibus analys of variance (ANOVA):\n"
		retstr += "Total squared deviation of subjects from the grand mean:%.3f\n"%(self.grandss)
		retstr += "Total squared deviation of group means from grand mean:%.3f\n"%(self.groupss)
		retstr += "Total squared deviation of subjects from their group mean:%.3f\n"%(self.subjectss)
		retstr += "Check that %.3f + %.3f = %.3f\n"%(self.groupss,self.subjectss,self.grandss)
		retstr += "%.3f + %.3f = %.3f\n"%(self.groupss,self.subjectss,self.groupss+self.subjectss)
		retstr += "Mean squared deviation (variance) of group means from grand mean (F Numerator):%.3f\n"%(self.groupms)
		retstr += "Mean squared deviation (variance) of subjects from their group mean (F Denominator):%.3f\n"%(self.subjectms)
		retstr += "Degrees of freedom numerator,denomenator:%s,%s\n"%(self.groupdf,self.subjectdf)
		retstr += "F=%.3f,p=%.5f"%(self.F,self.p)
		
		retstr += "\n\n"
		retstr += "Contrast analysis of variance (ANOVA):\n"
		for contrast_name in self.contrasts.keys():
			retstr += "Contrast:" + contrast_name + '\n'
			retstr += "Variance (1df so same as SS):%.3f\n"%(self.contrasts_ms[contrast_name])
			retstr += "Degrees of freedom numerator,denomenator:%3f,%3f\n"%(1,self.subjectdf)
			retstr += "F:%.3f,p=%.5f"%(self.contrasts_F[contrast_name],self.contrasts_p[contrast_name])
			retstr += "\n\n"
		
		
		
		return retstr