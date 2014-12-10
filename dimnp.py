#!/usr/bin/env python
#-*-coding:utf-8-*-

import os,sys
import time
import math

def rundimnp():
	if len(sys.argv) != 5:
		print 'Error'
		return 0
	ipathr= sys.argv[1]   #absolute
	opathr = sys.argv[2]   #absolute
	format = sys.argv[3]   #format:bed,bowtie
	pair_end = sys.argv[4]   # 0 False  1 True
	sigma = 35
	times_of_sigma = 2
	Nucl_signal,Nucl_peaksl = preprocess(ipathr,opathr,format,pair_end,sigma,times_of_sigma)


def preprocess(ipathr,opathr,format,pair_end,sigma,times_of_sigma):
	if not(os.path.exists(opathr)): os.mkdir(opathr)
	ipat = ipathr.split('*')
	Num = len(ipat)
	nucl_signal = {}
	nucl_peaksl = {}
	
	for ipath in ipat:
		tmp1 = os.path.split(ipath)
		tmp2 = os.path.splitext(tmp1[1])
		opath = os.path.join(opathr,tmp2[0])
		print ipath,opathr,opath
		if not(os.path.exists(opath)): os.mkdir(opath)
		nucl = Preprocess(ipath,opath,format,pair_end,sigma,times_of_sigma)
		name = ipath + '_signal'
		nucl_signal[name],nucl_peaksl[name]= nucl.runPreprocess()
	return nucl_signal,nucl_peaksl
	print nucl_signal,nucl_peaksl


class Preprocess:
	def __init__(self,ipath,opath,format,pair_end,sigma,times_of_sigma):
		self.ipath = ipath
		self.opath = opath
		self.format = format
		self.paired = pair_end
		self.sigma = sigma
		self.times_of_sigma = times_of_sigma
		self.chlength = {}
		self.nucl_profile = {}
		self.nucl_peaksl = {}
		self.nucl_peaksh = {}

	def make_profile(self):
		print 'Calculating Nuclesome Profile...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		opath1 = os.path.join(self.opath,'split')
		opath2 = os.path.join(self.opath,'profile')
		if not(os.path.exists(opath1)): os.mkdir(opath1)
		if not(os.path.exists(opath2)): os.mkdir(opath2)
		outfiles = {}  #step 1
		read_center = {}  #step 2   # include the reads in all chromesomes,seperate by key

		if self.format == 'bed':
			rawfile = open(self.ipath, 'r')
			idx = 1
			for line in rawfile:
				ch =  line.split()[0]
				outfile = outfiles.get(ch, None)
				if outfile is None:
					outfile = open(os.path.join(opath1, ch), 'w')
					outfiles[ch] = outfile
				outfile.write(line)
				# if ((idx % 1000000) == 0):
				# 	print "processed %i lines" % idx
				# idx += 1
			rawfile.close()
			for v in outfiles.values():
				v.close()
			for parent,dirnames,filenames in os.walk(opath1):  # get center
				for filename in filenames: 
					idx = 1
					if filename not in read_center:
						read_center[filename] = []
					if filename not in self.chlength:
						self.chlength[filename] = 0
					if filename not in self.nucl_profile:
						self.nucl_profile[filename] = []
					rawfile1 = open(os.path.join(opath1,filename),'r')
					for line in rawfile1:
						if line.split()[5] == '+': read_center[filename].append(int(line.split()[1]) + 75)
						elif line.split()[5] == '-': read_center[filename].append(int(line.split()[2]) - 75)
						# if ((idx % 100000) == 0): print "processed %i lines" % idx
						# idx += 1
					rawfile1.close()
			for key in read_center: # get length & initial self.nucl_profile
				self.chlength[key] = max(read_center[key]) + 100
				self.nucl_profile[key] = [0] * self.chlength[key]
				for i in read_center[key]:
					if i >= 37 and i <= self.chlength[key] - 37:
						idxs = xrange(i - 37, i + 37)
						for idx in idxs:
							self.nucl_profile[key][idx] += 1
			print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())				
			for key in self.nucl_profile:
				fr = open(os.path.join(opath2,key),'w')
				fr.writelines(str(self.nucl_profile[key]).strip('[]'))
				fr.close
			print 'self.chlength',self.chlength

		elif self.format == 'bowtie' and self.paired == '0': # for single data
			rawfile = open(self.ipath, 'r')
			idx = 1
			for line in rawfile:
				ch =  line.split()[4]
				outfile = outfiles.get(ch, None)
				if outfile is None:
					outfile = open(os.path.join(opath1, ch), 'w')
					outfiles[ch] = outfile
				outfile.write(line)
				# if ((idx % 100000) == 0):
				# 	print "processed %i lines" % idx
				# idx += 1
			rawfile.close()
			for v in outfiles.values():
				v.close()
			for parent,dirnames,filenames in os.walk(opath1):  # get center
				for filename in filenames: 
					idx = 1
					if filename not in read_center:
						read_center[filename] = []
					if filename not in self.chlength:
						self.chlength[filename] = 0
					if filename not in self.nucl_profile:
						self.nucl_profile[filename] = []
					rawfile1 = open(os.path.join(opath1,filename),'r')
					for line in rawfile1:
						if line.split()[3] == '+': read_center[filename].append(int(line.split()[5]) + 75)
						elif line.split()[3] == '-': read_center[filename].append(int(line.split()[5]) + int(line.split()[2][7:])- 75)
						# if ((idx % 100000) == 0): print "processed %i lines" % idx
						# idx += 1
					rawfile1.close()
			for key in read_center: # get length & initial self.nucl_profile
				self.chlength[key] = max(read_center[key]) + 100
				self.nucl_profile[key] = [0] * self.chlength[key]
				for i in read_center[key]:
					if i >= 37 and i <= self.chlength[key] - 37:
						idxs = xrange(i - 37, i + 37)
						for idx in idxs:
							self.nucl_profile[key][idx] += 1
			print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
			for key in self.nucl_profile:
				fr = open(os.path.join(opath2,key),'w')
				fr.writelines(str(self.nucl_profile[key]).strip('[]'))
				fr.close
			print 'self.chlength',self.chlength

		elif self.format == 'bowtie' and self.paired == '1': # for self.paired data
			rawfile = open(self.ipath,'r')
			idx = 1
			while True:
				line1 = rawfile.readline().split()
				line2 = rawfile.readline().split()
				if line1 == [] or line2 == []:
					break
				else:
					name1,stra1,ch1,start1,seqlen1 = line1[0],line1[2],line1[3],line1[4],len(line1[5])
					name2,stra2,ch2,start2,seqlen2 = line2[0],line2[2],line2[3],line2[4],len(line2[5])
					if name1 != name2:
						print 'pair error --- single end reads:\n',line1,'\n'
						line1 = rawfile.readline().split()
						# if ((idx % 100000) == 0): print "processed %i reads" % idx
						# idx += 1	
						continue
					else:
						if ch1 not in read_center:
							read_center[ch1] = []
						if stra1 != stra2:
							read_center[ch1].append((int(start1) + int(start2) + seqlen1)//2)
						if stra1 == stra2:
								print 'pair error --- reads from same strand:\n',line1,'\n',line2,'\n'
						# if ((idx % 100000) == 0): print "processed %i reads" % idx
						# idx += 1
			for key in read_center: # get length & initial self.nucl_profile
				self.chlength[key] = max(read_center[key]) + 100
				self.nucl_profile[key] = [0] * self.chlength[key]
				for i in read_center[key]:
					if i >= 37 and i <= self.chlength[key] - 37:
						idxs = xrange(i - 37, i + 37)
						for idx in idxs:
							self.nucl_profile[key][idx] += 1
			print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
			for key in self.nucl_profile:
				fr = open(os.path.join(opath2,key),'w')
				fr.writelines(str(self.nucl_profile[key]).strip('[]'))
				fr.close
			print 'self.chlength',self.chlength

	def rmclonal(self):   # step1:rmclonal   将冗余置为10*av（覆盖度）
		print 'Removing Clonal Reads...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		for key in self.nucl_profile:
			av = sum(self.nucl_profile[key])/float(self.chlength[key])
			print key,'av:',av 
			cut = av * 10
			for i in range(self.chlength[key]):
				if self.nucl_profile[key][i] > cut: self.nucl_profile[key][i] = cut
		opath3 = os.path.join(self.opath,'rmclonal')
		if not(os.path.exists(opath3)): os.mkdir(opath3)
		print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		for key in self.nucl_profile:
			fr = open(os.path.join(opath3,key),'w')
			fr.writelines(str(self.nucl_profile[key]).strip('[]'))
			fr.close

	def smooth(self):    # step2:smooth   对核小体信号做高斯卷积平滑
		print 'Data Smooth by Gaussian convolution......','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		def convolution(x,Gaussian,key):
			y=0
			if x>=self.sigma*self.times_of_sigma and x<=self.chlength[key]-self.sigma*self.times_of_sigma-1:
				for n in range(x-self.sigma*self.times_of_sigma,x+self.sigma*self.times_of_sigma+1):
					y=y+self.nucl_profile[key][n]*Gaussian[n-(x-self.sigma*self.times_of_sigma)]
			elif x<self.sigma*self.times_of_sigma:
				for n in range(0,x+self.sigma*self.times_of_sigma+1):
					y=y+self.nucl_profile[key][n]*Gaussian[self.sigma*self.times_of_sigma-x+n]
			elif x>self.chlength[key]-self.sigma*self.times_of_sigma-1:
				for n in range(x-self.sigma*self.times_of_sigma,self.chlength[key]):
					y=y+self.nucl_profile[key][n]*Gaussian[self.sigma*self.times_of_sigma-x+n]
			return y

		for key in self.nucl_profile:
			Gaussian=[]
			for x0 in range(-self.sigma*self.times_of_sigma,self.sigma*self.times_of_sigma+1):
				x=x0*(-1)
				Gaussian.append(math.exp(-x*x/(2.0*self.sigma*self.sigma))/math.sqrt(2*(math.pi)*self.sigma*self.sigma))
			for j in range(self.chlength[key]):
				self.nucl_profile[key][j] = convolution(j,Gaussian,key)
		print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		opath4 = os.path.join(self.opath,'smooth')
		if not(os.path.exists(opath4)): os.mkdir(opath4)
		for key in self.nucl_profile:
			fr = open(os.path.join(opath4,key),'w')
			fr.writelines(str(self.nucl_profile[key]).strip('[]'))
			fr.close

	def Fnor(self):   # step3:Fnor(Foldchange)   将核小体信号转为Foldchange（归一化）
		print 'Foldchange Normalization...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		for key in self.nucl_profile:
			av_ = sum(self.nucl_profile[key])/float(self.chlength[key])
			print key,'av_:',av_
			for i in range(self.chlength[key]):
				self.nucl_profile[key][i] = self.nucl_profile[key][i]/av_
		opath5 = os.path.join(self.opath,'normalization')
		if not(os.path.exists(opath5)): os.mkdir(opath5)
		print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		for key in self.nucl_profile:
			fr = open(os.path.join(opath5,key),'w')
			fr.writelines(str(self.nucl_profile[key]).strip('[]'))
			fr.close


	def Findpeaks(self):
		print 'Finding peaks...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		for key in self.nucl_profile:
			maxl = []
			maxh = []
			# first_der = [0]
			# for i in range(1,len(self.nucl_profile)-1):
			# 	first_der.append(self.nucl_profile[key][i+1] - self.nucl_profile[key][i-1])
			# first_der.append = [0]	
			for i in range(1,len(self.nucl_profile[key])-1):
				if self.nucl_profile[key][i] > self.nucl_profile[key][i+1] and self.nucl_profile[key][i] > self.nucl_profile[key][i-1]:
					maxl.append(i)
					maxh.append(self.nucl_profile[key][i])
			self.nucl_peaksl[key] = maxl
			self.nucl_peaksh[key] = maxh
		opath6 = os.path.join(self.opath,'findpeaks')
		if not(os.path.exists(opath6)): os.mkdir(opath6)
		print 'Writing into files...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
		for key in self.nucl_peaksl:
			fr = open(os.path.join(opath6,key),'w')
			fr.writelines(str(self.nucl_peaksl[key]).strip('[]'))
			fr.close


	def runPreprocess(self):
		self.make_profile()
		self.rmclonal()
		self.smooth()
		self.Fnor()
		self.Findpeaks()
		return self.nucl_profile,self.nucl_peaksl


rundimnp()
