#!/usr/bin/env python
#-*-coding:utf-8-*-

import os,sys
import time
from math import log10
from rpy2.robjects import r,FloatVector


# fn为列表，存放数据绝对路径，chr表示待动态识别的染色体号,opath为文件输出目录
def caldiff(fn = ['/vagrant/OUTPUT/013656/rm_sm_normalization/Scchr01','/vagrant/OUTPUT/013657/rm_sm_normalization/Scchr01',\
	'/vagrant/OUTPUT/013658/rm_sm_normalization/Scchr01'],chrNum = 'Scchr01',opath = './',\
	cutoff = 0.05, win_size = 36, step = 1):

	diffPoints = []
	P_values = [0]

	Num = len(fn)
	inpt = [FloatVector(r['scan'](file= name,what =r.double(0),sep = ","))*10 for name in fn]
	# print 'inpt:',inpt
	length = min([len(sample) for sample in inpt])
	print 'length:',length
	mean = [r.mean(sample)[0] for sample in inpt]
	print 'mean:',mean

	print 'calcularing diffPoints...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())

	for i in range(1,length-1,step):
		# print 'i:',i
		idx = 1
		wmean = get_wmean(win_size,length,i,inpt)
		p_value = calPvalue(mean,wmean,Num,method = 'chisq')
		P_values.append(-log10(p_value))
		if p_value <= cutoff: diffPoints.append(i)
		# if ((idx % 1000) == 0): print "processed %i points" % idx
		# 	idx += 1

	P_values.append(0)
	print 'len(P_values):',len(P_values)
	print 'finish cal_diff','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())


	fr = open(os.path.join(opath,'diffPoints_'+ chrNum),'w')
	fr.writelines(str(diffPoints).strip('[]'))
	fr.close()
	fr = open(os.path.join(opath,'P_values_'+ chrNum),'w')
	fr.writelines(str(P_values).strip('[]'))
	fr.close()

	# return P_values,diffPoints




def get_wmean(win_size,length,i,inpt):
	if i <= win_size/2: wmean = [r.mean(sample[0:2*i])[0] for sample in inpt]
	elif i > win_size/2 and i <= length-1-win_size/2: wmean = [r.mean(sample[i-win_size/2:i+win_size/2])[0] for sample in inpt]
	elif i > length-1-win_size/2 and i < length-1: wmean = [r.mean(sample[2*i-length+1:length-1])[0] for sample in inpt]
	# print 'wmean:',wmean
	return wmean


def calPvalue(mean,wmean,Num,method = 'chisq'): #mean = list,wmean = list
	data = r.matrix(FloatVector(mean + wmean),nrow = Num)
	if method == 'chisq':
		test = r['chisq.test'](data)
		p_value = test[test.names.index('p.value')][0]
		if p_value < 1e-10: p_value = 1e-10
	elif method == 'fisher':
		test = r['fisher.test'](data)
		p_value = test[test.names.index('p.value')][0]
		if p_value < 1e-10: p_value = 1e-10
	# print 'p_value:',p_value
	return p_value

