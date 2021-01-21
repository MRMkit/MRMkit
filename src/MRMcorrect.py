import collections
import os
import glob
import sys
import math
import statistics
from bisect import bisect_left
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.style as mplstyle
mplstyle.use([ 'fast'])
import numpy as np
print('NumPy',np.__version__)
print('Matplotlib',matplotlib.__version__)


import commonfn

param_set={
        "mzML_files",
        "batch_info",
        "batch_correction",
        "transition_list",
        "smoothing_sigma",
        }
param_dict=commonfn.read_param(param_set)

sig=float(param_dict.get("smoothing_sigma",10))

corr0=int(param_dict["batch_correction"]) #set drift correction?
transition=param_dict["transition_list"]
t_list=commonfn.read_trans_stat(transition)

ppfile='peak_picking.txt'
if os.path.isfile(ppfile):
    print(ppfile+' exists.')
    input('Press Enter/Return to resume')

miscdir='trans'

plot0=corr0 #plot 
starti=4
with open('quant_auc_rt.txt') as quant0:
    lsp1=quant0.readline().rstrip('\n').split('\t')
    blocklen=int((len(lsp1)-starti)/3)
    mzML_files=lsp1[starti:starti+blocklen]

batch_info_file=glob.glob(param_dict["batch_info"])[0]


batch_files=collections.defaultdict(list) #batch to files map
file_batch=dict()  #file to batch map
file_type=dict()  #file to type map
with open(batch_info_file) as bfile:
    next(bfile)
    for line in bfile:
        if line[0]!='#'and line.rstrip():
            lsp=line.rstrip('\n').split('\t')
            if not lsp[0].endswith('.mzML'):
                lsp[0]+='.mzML'
            batch_files[lsp[1]].append(lsp[0])
            file_batch[lsp[0]]=lsp[1]
            file_type[lsp[0]]=lsp[2].upper()

files_in_bf=set(vv for v in batch_files.values() for vv in v)
files_in_dir=set(mzML_files)

with open('missing_files.txt','w') as miss0:
    for f0 in mzML_files:
        if f0 not in files_in_bf:
            miss0.write("{} not in {}\n".format(f0,batch_info_file))
    for f0 in sorted(files_in_bf):
        if f0 not in files_in_dir:
            miss0.write("{} not in {}\n".format(f0,'directory'))

norm_dat=[]
unnorm_dat=[]
with open('quant_auc_rt.txt') as quant:
    lsp1=quant.readline().rstrip('\n').split('\t')
    file_ord=[]
    batch_l=[]
    type_l=[]
    for file0 in [vv for v in batch_files.values() for vv in v]:
        if file0 in files_in_dir:
            file_ord.append(mzML_files.index(file0))
            batch_l.append(file_batch[file0])
            type_l.append(file_type[file0])
    for line in quant:
        lsp=line.rstrip('\n').split('\t')
        unnorm_dat.append(lsp[:starti]+[float(lsp[starti+fo]) for fo in file_ord])
        if lsp[starti+blocklen]:
            norm_dat.append(lsp[:starti]+[float(lsp[starti+blocklen+fo]) for fo in file_ord])

def dnorm(x):
    return math.exp(-.5*(x/sig)**2)

largest_b=max(len(v) for v in batch_files.values())
dnorm_dict={ x:dnorm(x) for x in range(-largest_b,largest_b)}



def gp_reg(ss,st0,en0):
    xs=list(range(st0,en0))
    yy=[y for x,y in ss if np.isfinite(y)]
    if len(yy)<3:
        return xs,[0 for _ in xs]
    q1,q3=np.quantile(yy,[.1,.9])
    ub=q3+max(1,q3-q1)
    lb=q1-max(1,q3-q1)
    ss_=[(y,z) for y,z in ss if lb<z<ub]
    xx=[x for x,_ in ss_]
    yy=[x for _,x in ss_]

    if corr0:
        Efs_yy=[sum(dnorm_dict[xs1-xx2]*yy[xxj]for xxj,xx2 in enumerate(xx))/sum(dnorm_dict[xs1-xx2]for xx2 in xx) for xs1 in xs]
    else:
        Efs_yy=[0 for _ in xs]
    return xs,Efs_yy


colord={'BQC':'r','TQC':'b','SAMPLE':'k'}
colorl=[colord[x] for x in type_l]
change_pt=[ii+1 for ii in range(len(batch_l)-1)if batch_l[ii]!=batch_l[ii+1]]
adj_dat=[]

with PdfPages('run_seq_correction.pdf') as pdf0, PdfPages('run_seq_log2.pdf') as pdf1:
    seq0=list(range(len(unnorm_dat[0])-starti))
    seq0BQC=[x for x,y in zip(seq0,type_l)if y=='BQC']
    seq0TQC=[x for x,y in zip(seq0,type_l)if y=='TQC']
    seq0sam=[x for x,y in zip(seq0,type_l)if y=='SAMPLE']

    for nn,nd0 in enumerate(unnorm_dat[:]):
        print(nn,nd0[0])
        plt.figure(figsize=(18,2))
        plt.title(nd0[0]+'    log2-transformed peak area')
        for cp0 in change_pt:
            plt.axvline(x=cp0-.5,color='k',lw=1)
        nd_=nd0[starti:]
        log2nd=np.log2(nd_)
        log2nd_=[x for x,y in zip(log2nd,type_l)if y=='BQC']
        plt.scatter(seq0BQC,log2nd_,s=2,color='r',label='BQC')
        log2nd_=[x for x,y in zip(log2nd,type_l)if y=='TQC']
        plt.scatter(seq0TQC,log2nd_,s=2,color='b',label='TQC')
        log2nd_=[x for x,y in zip(log2nd,type_l)if y=='SAMPLE']
        plt.scatter(seq0sam,log2nd_,s=2,color='k',alpha=.2)

        plt.legend()
        q1,q3=np.quantile([x for x in log2nd if np.isfinite(x)],[.1,.9])
        ub=q3+max(1,q3-q1)
        lb=q1-max(1,q3-q1)
        plt.ylim(lb-.5,ub+.5)

        plt.tight_layout()
        pdf1.savefig()
        plt.close()

    seq0=list(range(len(norm_dat[0])-starti))
    seq0BQC=[x for x,y in zip(seq0,type_l)if y=='BQC']
    seq0TQC=[x for x,y in zip(seq0,type_l)if y=='TQC']
    seq0sam=[x for x,y in zip(seq0,type_l)if y=='SAMPLE']
    for nn,nd0 in enumerate(norm_dat[:]):
        print(nn,nd0[0])
        nd_=nd0[starti:]
        if plot0:
            fig, (ax0, ax1) = plt.subplots(2, sharex=True,figsize=(18, 4))
            ax0.set_title(nd0[0]+'    log2-transformed normalised peak area')
            ax1.set_title('time trend and batch correction')

            for cp0 in change_pt:
                ax0.axvline(x=cp0-.5,color='k',lw=1)
                ax1.axvline(x=cp0-.5,color='k',lw=1)

        adj0_l=[]
        for cp0,cp1 in zip([0]+change_pt,change_pt+[len(batch_l)]):
            seq1=seq0[cp0:cp1]
            nd1=np.log2(nd_[cp0:cp1])

            ss=[(y,z) for x,y,z in zip(type_l[cp0:cp1],seq1,nd1)if x!='TQC']
            xy=gp_reg(ss,cp0,cp1)

            if corr0 and plot0:
                ax0.plot(xy[0], xy[1], 'g-',lw=1)
            adj0_l.append([(pt-gk) for pt,gk in zip(nd1,xy[1])])

        if plot0:
            log2nd=np.log2(nd_)
            log2nd_=[x for x,y in zip(log2nd,type_l)if y=='BQC']
            ax0.scatter(seq0BQC,log2nd_,s=2,color='r',label='BQC')
            log2nd_=[x for x,y in zip(log2nd,type_l)if y=='TQC']
            ax0.scatter(seq0TQC,log2nd_,s=2,color='b',label='TQC')
            log2nd_=[x for x,y in zip(log2nd,type_l)if y=='SAMPLE']
            ax0.scatter(seq0sam,log2nd_,s=2,color='k',alpha=.2)
            ax0.legend()
            q1,q3=np.quantile([x for x in log2nd if np.isfinite(x)],[.1,.9])
            ub=q3+max(1,q3-q1)
            lb=q1-max(1,q3-q1)
            ax0.set_ylim(lb-.5,ub+.5)

        nd=[math.log2(x) for x,y in zip(nd_,type_l) if x>0 and y!='TQC']
        if corr0:
            allmed=(statistics.median(nd)if nd else 0)
        else:
            allmed=0

        adj_dat0=[]
        for adj0,cp0,cp1 in zip(adj0_l,[0]+change_pt,change_pt+[len(batch_l)]):
            seq1=seq0[cp0:cp1]
            adj0_=[x for x,y in zip(adj0,type_l[cp0:cp1]) if np.isfinite(x) and y!='TQC']

            if corr0:
                adjmed=(statistics.median(adj0_) if adj0_ else 0)
            else:
                adjmed=0

            nd1=[(x-adjmed+allmed if y!='TQC' else tqcy) for y,x,tqcy in zip(type_l[cp0:cp1],adj0,np.log2(nd_[cp0:cp1]))]

            adj_dat0.extend(nd1)
            if corr0 and plot0:
                ss=[(y,z) for x,y,z in zip(type_l[cp0:cp1],seq1,nd1)if x!='TQC']
                xy=gp_reg(ss,cp0,cp1)

                ax1.plot(xy[0], xy[1], 'g-',lw=1)

        if plot0:
            adj_dat0_=[x for x,y in zip(adj_dat0,type_l)if y=='BQC']
            ax1.scatter(seq0BQC,adj_dat0_,s=2,color='r')
            adj_dat0_=[x for x,y in zip(adj_dat0,type_l)if y=='TQC']
            ax1.scatter(seq0TQC,adj_dat0_,s=2,color='b')
            adj_dat0_=[x for x,y in zip(adj_dat0,type_l)if y=='SAMPLE']
            ax1.scatter(seq0sam,adj_dat0_,s=2,color='k',alpha=.2)
            q1,q3=np.quantile([x for x in log2nd if np.isfinite(x)],[.1,.9])
            ub=q3+max(1,q3-q1)
            lb=q1-max(1,q3-q1)
            ax1.set_ylim(lb-.5,ub+.5)

        adj_dat.append(np.power(2,adj_dat0))

        if plot0:
            plt.tight_layout()
            pdf0.savefig()
            plt.close()


filenames0=[mzML_files[fo] for fo in file_ord]
with open('batch_adjusted.txt','w') as adjf:
    adjf.write('filename\tbatch\ttype\t'+'\t'.join(tn for tn,*_ in norm_dat)+'\n')

    adjf.write('Q1\t\t\t'+'\t'.join(x[1] for x in norm_dat)+'\n')
    adjf.write('Q3\t\t\t'+'\t'.join(x[2] for x in norm_dat)+'\n')
    adjf.write('RT\t\t\t'+'\t'.join(x[3] for x in norm_dat)+'\n')

    for cn,f0 in enumerate(filenames0):
        adjf.write(f0)
        adjf.write('\t'+file_batch[f0])
        adjf.write('\t'+file_type[f0])
        adjf.write(''.join('\t'+str(adj0[cn]) for adj0 in adj_dat))
        adjf.write('\n')


starti=3

with open('batch_adjusted.txt') as my_op:
    mylsp1=my_op.readline().rstrip('\n').split('\t')
    blocklen=int((len(mylsp1)-starti))
    mylsp1=mylsp1[starti:starti+blocklen]
    mylspQ1=my_op.readline().rstrip('\n').split('\t')[starti:starti+blocklen]
    mylspQ3=my_op.readline().rstrip('\n').split('\t')[starti:starti+blocklen]
    mylspRT=my_op.readline().rstrip('\n').split('\t')[starti:starti+blocklen]
    my_BQCdict=collections.defaultdict(list)
    batchBQC=set()
    my_TQCdict=collections.defaultdict(list)
    batchTQC=set()
    my_SAMdict=collections.defaultdict(list)
    batchSAM=set()
    for line in my_op:
        lsp=line.rstrip('\n').split('\t')
        if lsp[2]=='BQC':
            batchBQC.add(lsp[1])
            for l,cpdname,Q1,Q3,RT in zip(lsp[starti:],mylsp1,mylspQ1,mylspQ3,mylspRT):
                my_BQCdict[lsp[1]+'\t'+cpdname+'\t'+Q1+'\t'+Q3].append(float(l))
        elif lsp[2]=='TQC':
            batchTQC.add(lsp[1])
            for l,cpdname,Q1,Q3,RT in zip(lsp[starti:],mylsp1,mylspQ1,mylspQ3,mylspRT):
                my_TQCdict[lsp[1]+'\t'+cpdname+'\t'+Q1+'\t'+Q3].append(float(l))
        elif lsp[2]=='SAMPLE':
            batchSAM.add(lsp[1])
            for l,cpdname,Q1,Q3,RT in zip(lsp[starti:],mylsp1,mylspQ1,mylspQ3,mylspRT):
                my_SAMdict[lsp[1]+'\t'+cpdname+'\t'+Q1+'\t'+Q3].append(float(l))
    for k,v in my_SAMdict.items():
        q1,q3=np.quantile(v,[.1,.9])
        ub=q3+max(1,q3-q1)
        lb=q1-max(1,q3-q1)
        my_SAMdict[k]=[vv for vv in v if lb<vv<ub]
    batchTQC=sorted(batchTQC)
    batchSAM=sorted(batchSAM & batchBQC)
    batchBQC=sorted(batchBQC)
    npeaks=collections.defaultdict(list)
    for cpdname,Q1,Q3,RT in zip(mylsp1,mylspQ1,mylspQ3,mylspRT):
        cpd,rt=cpdname.rsplit(' RT',1)
        npeaks[(cpd,Q1,Q3,RT)].append(rt)


CoVBQC=collections.defaultdict(list)
with open("BQC_CoV.txt",'w') as qctab:
    qctab.write('name\tQ1\tQ3\tRT\tave #peaks\t%CoV_all\t'+'\t'.join('%CoV_'+x for x in batchBQC)+'\n')
    for cpdname,Q1,Q3,RT in zip(mylsp1,mylspQ1,mylspQ3,mylspRT):
        qctab.write('\t'.join((cpdname,Q1,Q3,RT)))
        qctab.write('\t{}'.format(len(npeaks[(cpdname.rsplit(' RT',1)[0],Q1,Q3,RT)])))
        selquantBQC=[xx for b0 in batchBQC for xx in my_BQCdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]]
        cov_=100*np.std(selquantBQC)/np.mean(selquantBQC)
        qctab.write('\t{:.1f}'.format(cov_))
        CoVBQC[(cpdname.rsplit(' RT',1)[0],Q1,Q3)].append(format(cov_,'.1f'))
        for b0 in batchBQC:
            selquantBQC=my_BQCdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]
            qctab.write('\t{:.1f}'.format(100*np.std(selquantBQC)/np.mean(selquantBQC)))
        qctab.write('\n')

CoVTQC=collections.defaultdict(list)
with open("TQC_CoV.txt",'w') as qctab:
    qctab.write('name\tQ1\tQ3\tRT\tave #peaks\t%CoV_all\t'+'\t'.join('%CoV_'+x for x in batchTQC)+'\n')
    for cpdname,Q1,Q3,RT in zip(mylsp1,mylspQ1,mylspQ3,mylspRT):
        qctab.write('\t'.join((cpdname,Q1,Q3,RT)))
        qctab.write('\t{}'.format(len(npeaks[(cpdname.rsplit(' RT',1)[0],Q1,Q3,RT)])))
        selquantTQC=[xx for b0 in batchTQC for xx in my_TQCdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]]
        cov_=100*np.std(selquantTQC)/np.mean(selquantTQC)
        qctab.write('\t{:.1f}'.format(cov_))
        CoVTQC[(cpdname.rsplit(' RT',1)[0],Q1,Q3)].append(format(cov_,'.1f'))
        for b0 in batchTQC:
            selquantTQC=my_TQCdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]
            qctab.write('\t{:.1f}'.format(100*np.std(selquantTQC)/np.mean(selquantTQC)))
        qctab.write('\n')

Dratio=collections.defaultdict(list)
with open("D_ratio.txt",'w') as qctab:
    qctab.write('name\tQ1\tQ3\tRT\tave #peaks\t%D-ratio_all\t'+'\t'.join('%D-ratio_'+x for x in batchTQC)+'\n')
    for cpdname,Q1,Q3,RT in zip(mylsp1,mylspQ1,mylspQ3,mylspRT):
        qctab.write('\t'.join((cpdname,Q1,Q3,RT)))
        qctab.write('\t{}'.format(len(npeaks[(cpdname.rsplit(' RT',1)[0],Q1,Q3,RT)])))
        selquantSAM=[xx for b0 in batchSAM for xx in my_SAMdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]]
        selquantBQC=[xx for b0 in batchBQC for xx in my_BQCdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]]
        Dratio_=100*np.std(selquantBQC)/np.std(selquantSAM)
        qctab.write('\t{:.1f}'.format(Dratio_))
        Dratio[(cpdname.rsplit(' RT',1)[0],Q1,Q3)].append(format(Dratio_,'.1f'))
        for b0 in batchSAM:
            selquantSAM=my_SAMdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]
            selquantBQC=my_BQCdict[b0+'\t'+cpdname+'\t'+Q1+'\t'+Q3]
            qctab.write('\t{:.1f}'.format(100*np.std(selquantBQC)/np.std(selquantSAM)))
        qctab.write('\n')




t_name=[x.name for x in t_list]
class_rt=collections.defaultdict(list)
for (cpd,q1,q3,rt),rt_l in npeaks.items():
    if len(rt_l)==1:
        class_rt[t_list[bisect_left(t_name,cpd)].lclass].append((float(rt),float(rt_l[0])))


class_intercept=dict()
with PdfPages('expected_detected.pdf') as pdf0:
    fig, (ax0) = plt.subplots(1, sharex=True,figsize=(9, 9))

    for (size0,color0),(class0,rt0) in zip(itertools.product(range(9),'bgrcmyk'),class_rt.items()):
        xx=[x for x,_ in rt0]
        yy=[x for _,x in rt0]
        ax0.plot([min(xx+yy), max(xx+yy)], [min(xx+yy), max(xx+yy)],lw=1,ls=':',color='k')
        ax0.scatter(xx,yy,marker='${}$'.format(size0),color=color0,label=class0)
        ax0.set_xlabel('expected')
        ax0.set_ylabel('detected')
        ax0.legend()
    plt.tight_layout()
    pdf0.savefig()
    plt.close()
    ii=0
    for class0,rt0 in class_rt.items():
        if len(rt0)>=5:
            if ii%4==0:
                plt.figure(figsize=(9, 9))
            ax0=plt.subplot(2,2,ii%4+1)
            ax0.set_title(class0)
            xx=[x for x,_ in rt0]
            yy=[x for _,x in rt0]
            ax0.plot([min(xx+yy), max(xx+yy)], [min(xx+yy), max(xx+yy)],lw=1,ls=':',color='k')#,transform=ax0.transAxes
            ax0.scatter(xx,yy,s=5)
            xx=np.array(xx)
            
            ax0.plot(xx, np.median(yy-xx) + xx,lw=1, color='r', label='fitted line')
            class_intercept[class0]=np.median(yy-xx)

            ax0.legend()

            ax0.set_xlabel('expected RT')
            ax0.set_ylabel('detected RT')

            if ii%4==3:
                pdf0.savefig()
                plt.close()
            ii+=1
    if ii%4!=0:
        pdf0.savefig()
        plt.close()


if os.path.isfile(ppfile):
    print(ppfile+' exists, will not be overwriten')
    sys.exit()

with open(ppfile,'w') as usersRT:
    usersRT.write('class\tname\tQ1\tQ3\t%CoV(TQC)\t%CoV(BQC)\tD-ratio\texpectedRT\tdetectedRT\tusersRT\n')
    for (cpd,q1,q3,rt),npeak in sorted(npeaks.items()):
        with open(os.path.join(miscdir,commonfn.format_filename('trans2_{}Q1_{}Q3_{}.txt'.format(cpd,q1,q3)))) as trans2:
            next(trans2)
            rt_l=trans2.readline().split('\t')
            minRT,maxRT=float(rt_l[0]),float(rt_l[-1])
            qlen=(maxRT-minRT)/4
            q25,q75=minRT+qlen,maxRT-qlen

        class0=t_list[bisect_left([x.name for x in t_list],cpd)].lclass
        itc=class_intercept.get(class0,0)
        pred=min(npeak,key=lambda x:abs(itc+float(rt)-float(x)))
        if not q25<float(pred)<q75:
            pred=''
        usersRT.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(class0,cpd,q1,q3,', '.join(CoVTQC[(cpd,q1,q3)]),', '.join(CoVBQC[(cpd,q1,q3)]),', '.join(Dratio[(cpd,q1,q3)]),rt,', '.join(npeak),pred))
