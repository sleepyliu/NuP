import os,sys
import NuclPreprocess as nu
import CalDiffPoints as ca


def rundimnp(IPATH,OPATH,fmt,pair_end,cutoff,method):
    IPATH = IPATH.split(':')
    print 'IPATH:',IPATH,'\n','OPATH:',OPATH,'\n','format:',fmt,'\n','paired:',pair_end,'\n',\
        'cutoff:',cutoff,'\n','method:',method,'\n'
    print 'Start to identify nucleosome positioning dynamics...'

    FILE_NAMES = []

    # step 1:preprocess for each sample
    for ipath in IPATH:
        filename = os.path.splitext(os.path.split(ipath)[1])[0]
        opath = os.path.join(OPATH,filename+'.wig')
        FILE_NAMES.append(opath)
        print opath
        s = nu.NuclPreprocess(ipath,opath,fmt,pair_end)
        s.runPreprocess()

    # step 2:cal diff points
    opath = os.path.join(OPATH,'diffresult.wig')
    I = ca.calDiffPoints(FILE_NAMES,opath,cutoff,method)
    I.runcalDiffPoints()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Identification of nucleosome positioning dynamics in multiple samples.')
    parser.add_argument('IPATH',default=None,help="Paths to Nucleosome datas,seperated by ':'")
    parser.add_argument('-O','--OPATH',dest='opath',default='..',help='the output directory')
    parser.add_argument('-f','--format',dest='fmt',default='bowtie',choices=('bed', 'bowtie'),help='Only support bed or bowtie format for input files.')
    parser.add_argument('-p','--paired',dest='pair_end',default=1,type=int,help='set to 1 if the input data is paired-end reads.') 
    parser.add_argument('-c','--cutoff',dest='cutoff',default=0.01,type=float,help='the p-value cutoff for caldiff')
    parser.add_argument('-m','--method',dest='method',default='chisq.test',choices=('chisq.test',),help='the statistic method used for caldiff')
    args = parser.parse_args(sys.argv[1:])
    kwargs = vars(args)
    rundimnp(**kwargs)


