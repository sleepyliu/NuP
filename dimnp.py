#!/usr/bin/env python
#-*-coding:utf-8-*-

import os,sys
import preprocess as prep
import difftest as df


def rundimnp():
	if len(sys.argv) < 2:
		print 'Error : Please input parameters!'
		return 0
	ipathr = sys.argv[1]   # 要分析的样本数据(bed,bowtie或profile)路径，以*分隔
	opathr = sys.argv[2]   # 所有输出文件的根目录，包括样本间按染色体比较的差异结果以及每个样本的单独目录
	chrN = sys.argv[3]   # 'whole' or 'Scchr01' 要进行样本间差异分析的染色体号或全基因组,'whole'只有在isprep= 1(True) 时才有效
	cutoff_ = float(sys.argv[4])  # cutoff值 P-value
	isprep = sys.argv[5]   # 0 False  1 True 是否首先需要计算核小体profile

	win_size_ = 36
	step_ = 1

	if not(os.path.exists(opathr)): os.mkdir(opathr)

	if isprep == '1':
		format = sys.argv[6]   # format : bed,bowtie
		pair_end = sys.argv[7]   # 0 False  1 True
		keep = sys.argv[8]   # 0 False  1 True 是否保存中间文件，若为F，则只保存每条染色体最终profile和核小体位点文件，否则保存生成profile时预处理的各步文件
		if len(sys.argv) != 9:
			print 'Error : Wrong number of parameters!'
			return 0
		sigma = 35
		times_of_sigma = 2

		print 'ipath:',ipathr,'\n','opath:',opathr,'\n','format:',format,'\t','pair_end:',pair_end,'\t','keep:',keep,'\n',\
		'chrN:',chrN,'\t','cutoff:',cutoff_,'\t','win_size:',win_size_,'\t','step:',step_,'\n',\
		'sigma:',sigma,'\t','times_of_sigma:',times_of_sigma

		filenames,ipaths,opaths = get_fn_path(ipathr,opathr)  # 样本单独的文件名，输入路径和输出路径
		# print 'filenames:',filenames,'\n','ipaths:',ipaths,'\n','opaths:',opaths

		# 第一步预处理：
		# 将原始文件按染色体划分后生成每条染色体的profile，并对其进行去冗余，平滑，归一化处理得到最终的profile，
		# 然后对其进行峰值识别得到每个核小体中心坐标
		for i in range(len(filenames)):
			print 'Preprocessing sample',filenames[i],':'
			ipath = ipaths[i]
			opath = opaths[i]
			print 'ipath:',ipath,'\t','opath:',opath
			nucl = prep.Preprocess(ipath,opath,format,pair_end,keep,sigma,times_of_sigma)
			nucl.runPreprocess()

		# 第二步动态识别：
		chrNs = get_chrNs(chrN,opaths)
		for chrn in chrNs:
			print 'Identifying dynamic positions in',chrn,':'
			fns = [os.path.join(opath,'rm_sm_normalization',chrn) for opath in opaths]
			print 'fns:',fns
			df.caldiff(fn = fns,chrNum = chrn,opath = opathr,cutoff = cutoff_,win_size = win_size_,step = step_)

	elif isprep == '0':
		if len(sys.argv) != 6:
			print 'Error : Wrong number of parameters!'
			return 0

		print 'ipath:',ipathr,'\n','opath:',opathr,'\n',\
		'chrN:',chrN,'\t','cutoff:',cutoff_,'\t','win_size:',win_size_,'\t','step:',step_,'\n',\

		# 第二步动态识别：
		print 'Identifying dynamic positions in',chrN,':'
		fns = ipathr.split('*')
		print 'fns:',fns
		df.caldiff(fn = fns,chrNum = chrN,opath = opathr,cutoff = cutoff_,win_size = win_size_,step = step_)



def get_fn_path(ipathr,opathr):
	ipaths = ipathr.split('*')
	filenames = [os.path.splitext(os.path.split(ipath)[1])[0] for ipath in ipaths]
	opaths = [os.path.join(opathr,fn) for fn in filenames]
	return filenames,ipaths,opaths


def get_chrNs(chrN,opaths):
	if chrN == 'whole':
		chrNs = []
		path = os.path.join(opaths[0],'rm_sm_normalization')
		for parent,dirnames,filenames in os.walk(path):
			for filename in filenames:
				chrNs.append(filename)
	else:
		chrNs = [chrN]
	print 'chrNs:',chrNs
	return chrNs






rundimnp()

