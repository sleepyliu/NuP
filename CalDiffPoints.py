#!/usr/bin/env python

import os
import time
from math import log10
from rpy2.robjects import r,FloatVector

FILE_NAMES = ['/vagrant/OUTPUT/013656.wig',
'/vagrant/OUTPUT/013657.wig',
'/vagrant/OUTPUT/013658.wig']

# OPATH = '/vagrant/OUTPUT/Diff'
OPATH = '/vagrant/OUTPUT/diffresult.wig'


class calDiffPoints:
    def __init__(self,file_names=FILE_NAMES,opath=OPATH,cutoff = 0.01,method = 'chisq.test'):
        self.file_names = file_names
        self.opath = opath
        self.cutoff = cutoff
        # self.range = diffrange
        self.method = method

        self.signal = {}   # key = sample_name  value = {chr:profile}
        self.chrlist = []
        self.length = {}

        self.P_values = {}
        self.diffPoints = {}
        # self.diffZone = {}
        self.diffcenter = {}

    def loaddata(self):
        for filename in self.file_names:
            try:
                fr = open(filename)
                self.signal[filename] = {}
                for line in fr:
                    if 'chrom' in line:
                        ch = line[6:].strip('\n')
                        self.signal[filename].setdefault(ch,[])
                        self.chrlist.append(ch)
                    else:
                        self.signal[filename][ch].append(float(line.strip('\n')))
                fr.close()
                print 'load %s completed' %filename
            except IOError:
                print 'Cannot open file %s' %filename
                return 0

    def caldiffPvalues(self):
        print 'calcularing Statistic P_values(-log10)...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        for ch in self.chrlist:
            self.P_values.setdefault(ch,[])
            self.diffPoints.setdefault(ch,[])
            self.length[ch] = min([len(self.signal[filename][ch]) for filename in self.signal])
            mean = [float(sum(self.signal[filename][ch])/len(self.signal[filename][ch])) for filename in self.signal]
            for j in range(self.length[ch]):
                data = [self.signal[filename][ch][j] for filename in self.signal]
                data.extend(mean)
                test_data = r.matrix(FloatVector(data),ncol =2)
                test = r[self.method](test_data)
                if str(test).split()[-1] != 'NA':
                    p_value = float(str(test).split()[-1])
                    self.P_values[ch].append(-log10(p_value))
                    if p_value < self.cutoff: self.diffPoints[ch].append(j)
                else:
                    self.P_values[ch].append(0)
            print "%s caldiff finished" %ch 

    # def getDynamicZone(self):
    #     print 'get Dynamic Zone...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
    #     for ch in self.chrlist:
    #         self.diffZone.setdefault(ch,[])
    #         start = self.diffPoints[ch][0]
    #         for i in range(1,len(self.diffPoints[ch])):
    #             end = self.diffPoints[ch][i-1]
    #             if self.diffPoints[ch][i] - end <= 5:continue
    #             else:
    #                 if end - start >= self.range: 
    #                     for j in range(start,end):self.diffZone[ch].append(j)
    #                 else:
    #                     diffcenter = (start + end)/2
    #                     if diffcenter - self.range/2 >= 0 and diffcenter + self.range/2 <= self.length[ch]:
    #                         for j in range(diffcenter - self.range/2,diffcenter + self.range/2):self.diffZone[ch].append(j)
    #                 start = self.diffPoints[ch][i]

    def getDiffCenter(self):
        print 'get diff center...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        for ch in self.chrlist:
            self.diffcenter.setdefault(ch,[])
            start = self.diffPoints[ch][0]
            for i in range(1,len(self.diffPoints[ch])):
                end = self.diffPoints[ch][i-1]
                if self.diffPoints[ch][i] - end <= 5:continue
                else:
                    self.diffcenter[ch].append((start + end)/2)
                start = self.diffPoints[ch][i]


    def writefiles(self):
        # def save(name,file_name):
        #     try:
        #         opt = os.path.join(self.opath,name)
        #         fr = open(opt,'w')
        #         for ch in file_name:
        #             fr.write('chrom='+ch+'\n')
        #             for i in range(len(file_name[ch])):
        #                 fr.write(str(file_name[ch][i])+'\n')
        #         fr.close()
        #     except IOError:
        #         print 'Can not write into %s' %opt
        #         return 0

        # if not(os.path.exists(self.opath)): os.mkdir(self.opath)
        # save('P_values',self.P_values)
        # save('diffPoints',self.diffPoints)
        # save('diffZone',self.diffZone)

        print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        try:
            fr = open(self.opath,'w')
            fr.write('chrom'+'\t'+'position'+'\t'+'p_value(-log10)'+'\n')
            for ch in self.diffcenter:
                for i in range(len(self.diffcenter[ch])):
                    fr.write(ch+'\t'+str(self.diffcenter[ch][i])+'\t'+str(self.P_values[ch][i])+'\n')
            fr.close()
            print 'Completed.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        except IOError:
            print 'Can not write into %s' % self.opath
            return 0


    def runcalDiffPoints(self):
        self.loaddata()
        self.caldiffPvalues()
        # self.getDynamicZone()
        self.getDiffCenter()
        self.writefiles()

