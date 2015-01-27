#!/usr/bin/env python
#-*-coding:utf-8-*-

import os,sys
import preprocess as prep
import difftest as df


def rundimnp(ipathr, opathr, chrN, cutoff, isprep, winsize, step, format, pairend, keep):
    if not(os.path.exists(opathr)): os.mkdir(opathr)

    if isprep:
        sigma = 35
        times_of_sigma = 2

        print 'ipath:',ipathr,'\n','opath:',opathr,'\n','format:',format,'\t','pairend:',pairend,'\t','keep:',keep,'\n',\
        'chrN:',chrN,'\t','cutoff:',cutoff_,'\t','winsize:',winsize,'\t','step:',step,'\n',\
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
            nucl = prep.Preprocess(ipath,opath,format,pairend,keep,sigma,times_of_sigma)
            nucl.runPreprocess()

        # 第二步动态识别：
        chrNs = get_chrNs(chrN,opaths)
        for chrn in chrNs:
            print 'Identifying dynamic positions in',chrn,':'
            fns = [os.path.join(opath,'rm_sm_normalization',chrn) for opath in opaths]
            print 'fns:',fns
            df.caldiff(fn = fns,chrNum = chrn,opath = opathr,cutoff = cutoff_,winsize = winsize,step = step)

    elif isprep == '0':
        if len(sys.argv) != 6:
            print 'Error : Wrong number of parameters!'
            return 0

        print 'ipath:',ipathr,'\n','opath:',opathr,'\n',\
        'chrN:',chrN,'\t','cutoff:',cutoff_,'\t','winsize:',winsize,'\t','step:',step,'\n',\

        # 第二步动态识别：
        print 'Identifying dynamic positions in',chrN,':'
        fns = ipathr.split('*')
        print 'fns:',fns
        df.caldiff(fn = fns,chrNum = chrN,opath = opathr,cutoff = cutoff_,winsize = winsize,step = step)



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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Process.')
    parser.add_argument('-i', '--ipathr', nargs='+', help='要分析的样本数据(bed,bowtie或profile)路径，以空格分隔')
    parser.add_argument('-o', '--opathr', help='所有输出文件的根目录，包括样本间按染色体比较的差异结果以及每个样本的单独目录')
    parser.add_argument('-r', '--chrn', help="'whole' or 'Scchr01' 要进行样本间差异分析的染色体号或全基因组,'whole'只有在 isprep=True 时才有效")
    parser.add_argument('-c', '--cutoff', type=float, help='cutoff值 P-value')
    parser.add_argument('-p', '--isprep', action='store_true', help='是否首先需要计算核小体profile')
    parser.add_argument('-s', '--winsize', type=int, default=36, help='窗口大小')
    parser.add_argument('-t', '--step', type=int, default=1, help='步长')
    parser.add_argument('-f', '--format', choices=('bed', 'bowtie'), help='格式，bed or format')
    parser.add_argument('-e', '--pairend', action='store_true', help='是否有pair end')
    parser.add_argument('-k', '--keep', action='store_true', help='是否保存中间文件，若为False，则只保存每条染色体最终profile和核小体位点文件，否则保存生成profile时预处理的各步文件')

    args = parser.parse_args(sys.argv[1:])
    kwargs = vars(args)
    rundimnp(**kwargs)
