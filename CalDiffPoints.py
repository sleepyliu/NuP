#!/usr/bin/env python

import os
import time
from math import log10
from rpy2.robjects import r,FloatVector

FILE_NAMES = ['/vagrant/OUTPUT/013656/Scchr04',
'/vagrant/OUTPUT/013657/Scchr04',
'/vagrant/OUTPUT/013658/Scchr04']

OPATH = '/vagrant/OUTPUT/Diff/Scchr04'


class calDiffPoints:
    def __init__(self,file_names=FILE_NAMES,opath=OPATH,cutoff = 0.01,diffrange = 150,method = 'chisq.test'):
        self.file_names = file_names
        self.opath = opath
        self.cutoff = cutoff
        self.range = diffrange
        self.method = method

        self.signal = {}
        self.length = 0
        self.P_values = []
        self.diffPoints = []
        self.diffZone = []

    def loaddata(self):
        for i in self.file_names:
            try:
                print 'load %s' %i
                self.signal[i] = r['scan'](file=i,what =r.double(0),sep = ",")
            except IOError:
                print 'Cannot open file %s' %i
                return 0

    def caldiffPvalues(self):
        self.length = min([len(self.signal[i]) for i in self.signal])
        mean = [float(sum(self.signal[i])/len(self.signal[i])) for i in self.signal]

        print 'calcularing Statistic P_values(-log10)...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        for j in range(self.length):
            data = [self.signal[i][j] for i in self.signal]
            data.extend(mean)
            test_data = r.matrix(FloatVector(data),ncol =2)
            test = r[self.method](test_data)
            if str(test).split()[-1] != 'NA':
                p_value = float(str(test).split()[-1]) 
                self.P_values.append(-log10(p_value))
                if p_value < self.cutoff: self.diffPoints.append(j)
            else:
                self.P_values.append(0)
            if (j%1000 ==0):print "processed %i points" % j

    def getDynamicZone(self):
        start = self.diffPoints[0]
        for i in range(1,len(self.diffPoints)):
            end = self.diffPoints[i-1]
            if self.diffPoints[i] - end <= 5:continue
            else:
                if end - start >= self.range: 
                    for j in range(start,end):self.diffZone.append(j)
                else:
                    diffcenter = (start + end)/2
                    if diffcenter - self.range/2 >= 0 and diffcenter + self.range/2 <= self.length:
                        for j in range(diffcenter - self.range/2,diffcenter + self.range/2):self.diffZone.append(j)
                start = self.diffPoints[i]

    def writefiles(self):
        if not(os.path.exists(self.opath)): os.mkdir(self.opath)
        try:
            fr = open(os.path.join(self.opath,'P_values'),'w')
            fr.writelines(str(self.P_values).strip('[]'))
            fr.close()
            fr = open(os.path.join(self.opath,'diffPoints'),'w')
            fr.writelines(str(self.diffPoints).strip('[]'))
            fr.close()            
            fr = open(os.path.join(self.opath,'diffZone'),'w')
            fr.writelines(str(self.diffZone).strip('[]'))
            fr.close()
        except IOError:
            print 'Cannot write into file %s' %self.opath
            return 0

    def runcalDiffPoints(self):
        self.loaddata()
        self.caldiffPvalues()
        self.getDynamicZone()
        self.writefiles()

