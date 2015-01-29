#!/usr/bin/env python

import os
import time
import math


class nuclPreprocess:
    def __init__(self,ipath,opath,fmt,pair_end):
        self.ipath = ipath
        self.opath = opath
        self.fmt = fmt
        self.paired = pair_end

        self.chlength = {}
        self.nucl_profile = {}

    def make_profile(self):  # step1
        print 'Calculating Nuclesome Profile from %s ...' % self.ipath,'\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        try:
            rawfile = open(self.ipath, 'r')
        except IOError:
            print 'Can not open file %s' % self.ipath
            return 0

        read_center = {}  # key = chr,value = center of each read

        if self.fmt == 'bed':
            for line in rawfile:
                ch,start,dire =  line.split()[0],line.split()[1],line.split()[5]
                read_center.setdefault(ch,[])
                if dire == '+':read_center[ch].append(int(start) + 75)
                elif dire == '-':read_center[ch].append(int(start) - 75)

        elif self.fmt == 'bowtie' and self.paired == 0:
            for line in rawfile:
                ch,start,dire =  line.split()[4],line.split()[5],line.split()[3]
                read_center.setdefault(ch,[])
                if dire == '+':read_center[ch].append(int(start) + 75)
                elif dire == '-':read_center[ch].append(int(start) - 75)

        elif self.fmt == 'bowtie' and self.paired == 1:
            lname, ldire, lstart = '','',''
            for line in rawfile:
                splited = line.split()
                name, dire, ch, start, seqlen = splited[0],splited[2],splited[3],splited[4],len(splited[5])
                if name == lname and dire != ldire:
                    read_center.setdefault(ch, [])
                    read_center[ch].append((int(start) + int(lstart) + seqlen)//2)
                lname, ldire, lstart = name, dire, start

        rawfile.close()

        # calculate length and initial self.nucl_profile
        for key in read_center:
            self.chlength[key] = max(read_center[key]) + 100
            self.nucl_profile[key] = [0] * self.chlength[key]
            for i in read_center[key]:
                if i >= 37 and i <= self.chlength[key] - 37:
                    idxs = xrange(i - 37, i + 37)
                    for idx in idxs:
                        self.nucl_profile[key][idx] += 1
        # print 'self.chlength',self.chlength


    def rmclonal(self):  # step2:rmclonal(>10*av)
        print 'Removing Clonal Reads...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        for key in self.nucl_profile:
            av = sum(self.nucl_profile[key])/float(self.chlength[key])
            cut = av * 10
            for i in range(self.chlength[key]):
                if self.nucl_profile[key][i] > cut: self.nucl_profile[key][i] = cut


    def smooth(self):    # step3:Gaussian smooth
        print 'Data Smooth by Gaussian convolution...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        sigma = 35
        times_of_sigma = 2
        def convolution(x,Gaussian,key):
            y=0
            if x>=sigma*times_of_sigma and x<=self.chlength[key]-sigma*times_of_sigma-1:
                for n in range(x-sigma*times_of_sigma,x+sigma*times_of_sigma+1):
                    y=y+self.nucl_profile[key][n]*Gaussian[n-(x-sigma*times_of_sigma)]
            elif x<sigma*times_of_sigma:
                for n in range(0,x+sigma*times_of_sigma+1):
                    y=y+self.nucl_profile[key][n]*Gaussian[sigma*times_of_sigma-x+n]
            elif x>self.chlength[key]-sigma*times_of_sigma-1:
                for n in range(x-sigma*times_of_sigma,self.chlength[key]):
                    y=y+self.nucl_profile[key][n]*Gaussian[sigma*times_of_sigma-x+n]
            return y

        for key in self.nucl_profile:
            Gaussian=[]
            for x0 in range(-sigma*times_of_sigma,sigma*times_of_sigma+1):
                x=x0*(-1)
                Gaussian.append(math.exp(-x*x/(2.0*sigma*sigma))/math.sqrt(2*(math.pi)*sigma*sigma))
            for j in range(self.chlength[key]):
                self.nucl_profile[key][j] = convolution(j,Gaussian,key)


    def Fnor(self):   # step4:Fnor(Foldchange)
        print 'Foldchange Normalization...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        for key in self.nucl_profile:
            av_ = sum(self.nucl_profile[key])/float(self.chlength[key])
            # print key,'av_:',av_
            for i in range(self.chlength[key]):
                self.nucl_profile[key][i] = self.nucl_profile[key][i]/av_


    def set_normalization_level(self):
        av_set = 30
        for key in self.nucl_profile:
            for i in range(self.chlength[key]):
                self.nucl_profile[key][i] = self.nucl_profile[key][i] * av_set


    def writefiles(self):
        print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        try:
            fr = open(self.opath,'w')
            for key in self.nucl_profile:
                fr.write('chrom='+key+'\n')
                for i in range(len(self.nucl_profile[key])):
                    fr.write(str(self.nucl_profile[key][i])+'\n')
            fr.close()
            print 'Completed.','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
        except IOError:
            print 'Can not write into %s' % self.opath
            return 0


    def runPreprocess(self):
        self.make_profile()
        self.rmclonal()
        self.smooth()
        self.Fnor()
        self.set_normalization_level()
        self.writefiles()
        # return self.nucl_profile
        
