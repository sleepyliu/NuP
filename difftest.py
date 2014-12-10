#!/usr/bin/env python

import os,sys
import time
from math import log10
from rpy2.robjects import r

in1 = r['scan'](file='/vagrant/OUTPUT/output_013656/normalization/Scchr01',what =r.double(0),sep = ",")
in2 = r['scan'](file='/vagrant/OUTPUT/output_013657/normalization/Scchr01',what =r.double(0),sep = ",")
in3 = r['scan'](file='/vagrant/OUTPUT/output_013658/normalization/Scchr01',what =r.double(0),sep = ",")

cutoff = 0.1
win_size = 36
step = 1
diffPoints = []
P_values_chisq = [1]
P_values_fisher = [1]

length = min(len(in1),len(in2),len(in3))
r1 = r.mean(in1)
r2 = r.mean(in2)
r3 = r.mean(in3)
print r1,r2,r3
print 'length:',length

print 'calcularing diffPoints...','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())
idx = 2
i = 1
while i<= win_size/2:
	i1 = r.mean(in1[0:2*i])
	i2 = r.mean(in2[0:2*i])
	i3 = r.mean(in3[0:2*i])
	data = r.matrix(r.c(i1,r1,i2,r2,i3,r3),nrow = 2)
	test = r['chisq.test'](data)
	p_value = float(str(test).split()[21])
	P_values_chisq.append(-log10(p_value))
	test2 = r['fisher.test'](data)
	p_value2 = float(str(test2).split()[9])
	P_values_fisher.append(-log10(p_value2))
	# if p_value < cutoff: diffPoints.append(i)
	if ((idx % 1000) == 0): print "processed %i points" % idx
	idx += 1
	i += step

while i > win_size/2 and i <= length-1-win_size/2:
	i1 = r.mean(in1[i-win_size/2:i+win_size/2])
	i2 = r.mean(in2[i-win_size/2:i+win_size/2])
	i3 = r.mean(in3[i-win_size/2:i+win_size/2])
	data = r.matrix(r.c(i1,r1,i2,r2,i3,r3),nrow = 2)
	test = r['chisq.test'](data)
	p_value = float(str(test).split()[21])
	P_values_chisq.append(-log10(p_value))
	test2 = r['fisher.test'](data)
	p_value2 = float(str(test2).split()[9])
	P_values_fisher.append(-log10(p_value2))
	# if p_value < cutoff: diffPoints.append(i)
	if ((idx % 1000) == 0): print "processed %i points" % idx
	idx += 1
	i += step

while i > length-1-win_size/2 and i < length-1:
	i1 = r.mean(in1[2*i-length+1:length-1])
	i2 = r.mean(in2[2*i-length+1:length-1])
	i3 = r.mean(in3[2*i-length+1:length-1])
	data = r.matrix(r.c(i1,r1,i2,r2,i3,r3),nrow = 2)
	test = r['chisq.test'](data)
	p_value = float(str(test).split()[21])
	P_values_chisq.append(-log10(p_value))
	test2 = r['fisher.test'](data)
	p_value2 = float(str(test2).split()[9])
	P_values_fisher.append(-log10(p_value2))
	# if p_value < cutoff: diffPoints.append(i)
	if ((idx % 1000) == 0): print "processed %i windows" % idx
	idx += 1       
	i += step
	
P_values_chisq.append(0)
P_values_fisher.append(0)
print 'len(P_values_chisq):',len(P_values_chisq)
print 'len(P_values_fisher):',len(P_values_fisher)
print 'finish cal_diff','\t',time.strftime('%Y-%m-%d %A %H:%M:%S', time.localtime())

# output = './diffPoints'
# fr = open(output,'w')
# fr.writelines(str(diffPoints).strip('[]'))
# fr.close()

output = './P_values_chisq'
fr = open(output,'w')
fr.writelines(str(P_values_chisq).strip('[]'))
fr.close()

output = './P_values_fisher'
fr = open(output,'w')
fr.writelines(str(P_values_fisher).strip('[]'))
fr.close()
