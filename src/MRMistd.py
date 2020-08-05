import operator
import math
from bisect import bisect_left
import sys
import concurrent.futures
import os
import glob
import itertools
import time
from functools import partial


import commonfn
import MRMcwt
import MRMeic


param_set={
        "mzML_files",
        "transition_list",
        "length_of_ion_chromatogram",
        "RT_difference",
        }
param_dict=commonfn.read_param(param_set)

transition=param_dict["transition_list"]

RTdiff=param_dict.get("RT_difference")
RTdiff=float(RTdiff) if RTdiff is not None else None
if RTdiff is not None:
    print('RTdiff',RTdiff)

t_list=commonfn.read_trans_stat(transition)
istd_dict=commonfn.get_istd_dict(t_list)
print(len(istd_dict),'ISTDs')
t_list=[x for x in t_list if x.istd in istd_dict]

mzML_files=[]
for x in param_dict["mzML_files"].split():
    mzML_files.extend(glob.glob(x))
mzML_files.sort()
mzML_files=mzML_files[:]
print(len(mzML_files),'mzML files')

miscdir='trans'
if os.path.exists(miscdir):
    print('rerun?')
    input('Press Enter/Return to continue')
else:
    os.mkdir(miscdir)
trans2filename={trans:os.path.join(miscdir,commonfn.format_filename('trans_{}Q1_{:.1f}Q3_{:.1f}.txt'.format(trans.name,trans.q1,trans.q3))) for trans in t_list}

peak_w=[int(x) for x in param_dict["length_of_ion_chromatogram"].split()]
wave_scales=[x/2 for x in list(range(peak_w[1]+2,peak_w[0]-3,-2))]
wave_scales=[2,2.5]+list(range(3,20))+[20*1.05**i for i in range(20)]
ws0=bisect_left(wave_scales,peak_w[0]/2)
ws1=bisect_left(wave_scales,peak_w[1]/2)

wave_scales=wave_scales[:ws1+1+1]#[::-1]
wave_sqrt=[math.sqrt(x) for x in wave_scales]
findridge=partial(MRMcwt.findridge0,wave_scales,wave_sqrt,False)
findridge2=partial(MRMcwt.findridge0,wave_scales,wave_sqrt,True)


for trans in t_list:
    open(trans2filename[trans],'w').write('{}\n'.format(trans.name))

def write_to_trans(eic_dict,trans):

    rt_f=commonfn.find_trans([*eic_dict],trans)
    if rt_f is None:
        return
    rt_I=eic_dict[rt_f]
    if RTdiff is not None:# and rt_I[0][-1]-rt_I[0][0]>200:
        rt_I=[(x,y) for x,y in zip(*rt_I) if abs(trans.rt-x)<RTdiff]
        rt_I=[x for x,_ in rt_I],[x for _,x in rt_I]

    if trans.name == trans.istd:
        rt_,I_=rt_I
        firstdecile=sorted(I_)[int(len(I_)/10)]
        I_=[max(0,i-firstdecile) for i in I_]

        rerunw=3*(rt_f.maxRT-rt_f.minRT)/(len(rt_)-1)
        feats=findridge(rerunw,list(zip(rt_,I_)))
        feats_=findridge2(rerunw,list(zip(rt_,I_)))
        feats_=[p for p in feats_ if feats[0].begin<p.rt<feats[0].end]
        if len(feats_)>1:
            rt_l=sorted(x.rt for x in feats_)
            postop=bisect_left(rt_l,feats_[0].rt)
            if 0<postop<len(rt_l)-1:
                begin0=(rt_l[postop-1]+rt_l[postop])/2
                end0=(rt_l[postop]+rt_l[postop+1])/2
            elif postop==0:
                begin0=feats[0].begin
                end0=(rt_l[postop]+rt_l[postop+1])/2
            elif postop==len(rt_l)-1:
                begin0=(rt_l[postop-1]+rt_l[postop])/2
                end0=feats[0].end
            pos0=bisect_left(rt_,begin0)
            pos1=bisect_left(rt_,end0)
            RT_I_=list(zip(rt_,I_))[pos0:pos1]
            auc=0
            feats=[commonfn.Peak(feats_[0].rt,feats_[0].sc*2,None,None,auc,begin0,end0)]
    else:
        feats=[]
    return rt_I,feats



start_time = time.time()


with open('trans_mzML_list.txt','w')as mf_list:
    for mf in mzML_files:
        basename0=os.path.basename(mf)
        mf_list.write(basename0+'\n')

open('missing_compounds.txt','w')

nfile=math.ceil(len(mzML_files)/math.ceil(len(mzML_files)/300))
for sta_ in range(0,len(mzML_files),nfile):
    eic_dict_all=list(map(MRMeic.print_eic,mzML_files[sta_:sta_+nfile]))
    
    print('{} - {}'.format(sta_+1,sta_+nfile))
    def write_block(trans):
        with open(trans2filename[trans],'a') as fo,open('missing_compounds.txt','a')as miss_cpd:
            for mzml,eic_dict in zip(mzML_files[sta_:sta_+nfile],eic_dict_all):
                rt_I_feats=write_to_trans(eic_dict,trans)
                if rt_I_feats is None:
                    miss_cpd.write('"{}" not in "{}"\n'.format(trans.name,mzml))
                    continue
                rt_I,feats=rt_I_feats
                fo.write('{}\n'.format('\t'.join(str(x) for x in rt_I[0])))
                fo.write('{}\n'.format('\t'.join(str(x) for x in rt_I[1])))
                for feat in feats:
                    fo.write('{:.3f}\t{}\t{:.3f}\t{:.3f}\t{:.3f}\n'.format(feat.rt,feat.auc,feat.begin,feat.end,feat.sc))
                fo.write('\n')
    list(map(write_block,t_list))
    
miss_cpd=open('missing_compounds.txt','a')
for trans_file in glob.glob(os.path.join(miscdir,'trans_*Q1_*Q3_*.txt')):
    with open(trans_file) as trans_f:
        cs=0
        for line in trans_f:
            if line=='\n':
                cs+=1
        trans_f.seek(0)
        trans_name=trans_f.readline().rstrip()
    if cs<len(mzML_files):
        print('removed {}. Present in {}/{} samples'.format(trans_name,cs,len(mzML_files)))
        miss_cpd.write('removed {}. Present in {}/{} samples\n'.format(trans_name,cs,len(mzML_files)))
        os.remove(trans_file)


print("Run time = {:.1f} mins  ".format(((time.time() - start_time)/60)))
