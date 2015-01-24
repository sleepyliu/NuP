#!/usr/bin/env python

import time
from math import log10
from rpy2.robjects import r,FloatVector

FILE_NAMES = ('/vagrant/OUTPUT/013656/rm_sm_normalization/Scchr01',
'/vagrant/OUTPUT/013657/rm_sm_normalization/Scchr01',
'/vagrant/OUTPUT/013658/rm_sm_normalization/Scchr01')

OPATH = '/vagrant/new/Scchr01_diff'


class calDiffPoints:
    def __init__(self,file_names=FILE_NAMES,opath=OPATH,cutoff = 0.1,method = 'chisq.test'):
        self.file_names = file_names
        self.opath = opath
        self.cutoff = cutoff
        self.method = method

        self.signal = {}
        self.P_values = []
        self.diffPoints = []

    def loaddata(self):
        for i in self.file_names:
            try:
                print 'load %s' %i
                self.signal[i] = r['scan'](file=i,what =r.double(0),sep = ",")
            except IOError:
                print 'Cannot open file %s' %i
                return 0

    def caldiff(self):
        length = min([len(self.signal[i]) for i in self.signal])
        mean = [float(sum(self.signal[i])/len(self.signal[i])) for i in self.signal]

        print 'calcularing diffPoints...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        for j in range(length):
            data = [self.signal[i][j] for i in self.signal]
            data.extend(mean)
            test_data = r.matrix(FloatVector(data),ncol =2)
            test = r[self.method](test_data)  # calculate p_value
            p_value = float(str(test).split()[22])
            self.P_values.append(-log10(p_value))
            if p_value < self.cutoff: self.diffPoints.append(j)
            if (j%1000 ==0):print "processed %i points" % j

    def writefiles(self):
        try:
            fr = open(self.opath,'w')
            fr.writelines(str(self.P_values).strip('[]'))
            fr.close()
        except IOError:
            print 'Cannot write into file %s' %self.opath
            return 0

    def runcalDiffPoints(self):
        self.loaddata()
        self.caldiff()
        self.writefiles()
