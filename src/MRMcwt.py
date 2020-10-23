import collections
import operator
import math
from bisect import bisect_left
import sys
import concurrent.futures
import os
import glob
from multiprocessing import freeze_support
import time

import commonfn

Peak=commonfn.Peak


param_dict=commonfn.read_param({ "length_of_ion_chromatogram" })

peak_w=[int(x) for x in param_dict["length_of_ion_chromatogram"].split()]


Coef=collections.namedtuple('Coef',('rt sc coef'))


def findridge0(wave_scales,wave_sqrt,rerun,rerunw,EIC):
    eic_dict={rt:i for rt,i in EIC}
    rt_all=sorted(eic_dict)
    I_sub=sorted(eic_dict.values(),reverse=True)
    I_cut=I_sub[int(len(I_sub)*.6)]
    eic_rt=set()
    for x,y in eic_dict.items():
        pos=bisect_left(rt_all,x)
        if rt_all[0]<x<rt_all[-1] and y>I_cut:
            eic_rt.add((rt_all[pos-1]+x)/2)
            eic_rt.add((rt_all[pos+1]+x)/2)
    eic_rt=sorted(eic_rt)
    

    coefs = [[0]*len(eic_rt) for i in wave_scales]
    for xx,wave_scale in enumerate(wave_scales):
        bd=rerunw if rerun else max(wave_scale,rerunw)
        for yy,wave_loc in enumerate(eic_rt):
            pos0=bisect_left(rt_all,wave_loc-bd)
            pos1=bisect_left(rt_all,wave_loc+bd)
            rt_=rt_all[pos0:pos1]
            int_I=[0]*len(rt_)
            for i,rt0 in enumerate(rt_):
                if eic_dict[rt0]>0:# and rt0!=wave_loc:
                    tsig2=((rt0 - wave_loc)/wave_scale)**2
                    int_I[i]=eic_dict[rt0]*math.exp(-tsig2/2)*(1.-tsig2)
            coefs[xx][yy]=sum((I0+I1)*(rt1-rt0) for rt0,rt1,I0,I1 in zip(rt_,rt_[1:],int_I,int_I[1:]))/2/wave_sqrt[xx]


    rtdict={rt_i:n for n,rt_i in enumerate(eic_rt)}
    max_map = [[0]*len(eic_rt) for i in wave_scales]
    for xx,wave_scale in enumerate(wave_scales):
        max_coefs=[]
        coef_xx=[(eic_rt[n],x) for n,x in enumerate(coefs[xx])]
        while any(y>0 for _,y in coef_xx):
            max_coef=max(coef_xx,key=operator.itemgetter(1))
            max_i=coef_xx.index(max_coef)
            if max_i==0 or max_i==len(coef_xx)-1 or coef_xx[max_i-1][1]==0 or coef_xx[max_i+1][1]==0:
                coef_xx[max_i]=(coef_xx[max_i][0],0)
                continue
            max_coefs.append(max_coef)
            for i in range(len(coef_xx)):
                if abs(max_coef[0]-coef_xx[i][0])<(rerunw if rerun else wave_scale):
                    coef_xx[i]=(coef_xx[i][0],0)
        for rt,coef in max_coefs:
            max_map[xx][rtdict[rt]]=coef


    ridgelines=[]
    for xx,scale_coef in enumerate(max_map):
        for yy,coef in enumerate(scale_coef):
            added=False
            if coef:
                for rl in ridgelines:
                    if ((abs(bisect_left(rt_all,rl[-1].rt)-bisect_left(rt_all,eic_rt[yy]))<2 or abs(rl[-1].rt-eic_rt[yy])<1) and \
                            rl[-1].sc==wave_scales[xx-1]):
                        rl.append(Coef(eic_rt[yy],wave_scales[xx],coef))
                        added=True
                        break
                if not added:
                    ridgelines.append([Coef(eic_rt[yy],wave_scales[xx],coef)])

    ridgelines=[x for x in ridgelines if len(x)>4]

    if not rerun:
        for rd in ridgelines[:]:
            peak_loc=max(rd,key=operator.attrgetter('coef'))
            if not peak_w[0]<=peak_loc.sc*2<=peak_w[1]:
                ridgelines.remove(rd)
                continue

            pos0=bisect_left(rt_all,peak_loc.rt-peak_loc.sc)
            pos1=bisect_left(rt_all,peak_loc.rt+peak_loc.sc)
            if pos0==0 or pos1==len(rt_all):# or sum(1 for rt0 in rt_all[pos0:pos1] if rt0 in eic_dict)/(pos1-pos0)<.5:
                ridgelines.remove(rd)
                continue

    peaks=[]


    if not rerun:
        for rd in ridgelines:
            peak_loc=max(rd,key=operator.attrgetter('coef'))

            if peak_loc.sc>len(rd) and rd.index(peak_loc)<=len(rd)-5:
                continue

            pos0=bisect_left(rt_all,peak_loc.rt-peak_loc.sc-1)
            pos1=bisect_left(rt_all,peak_loc.rt+peak_loc.sc+1)

            posm=bisect_left(rt_all,peak_loc.rt)
            pos2=bisect_left(rt_all,peak_loc.rt-2.5*peak_loc.sc)
            pos3=bisect_left(rt_all,peak_loc.rt+2.5*peak_loc.sc)
            apex_h=(eic_dict[rt_all[posm-1]]+eic_dict[rt_all[posm]])/12
            if pos2<pos0-1 and max(eic_dict[rt0]for rt0 in rt_all[pos2:pos0-1])<apex_h:
                pos0=pos2
            if pos1+1<pos3 and max(eic_dict[rt0]for rt0 in rt_all[pos1+1:pos3])<apex_h:
                pos1=pos3
            auc=sum((eic_dict[rt0]+eic_dict[rt1])*(rt1-rt0) for rt0,rt1 in zip(rt_all[pos0:],rt_all[pos0+1:pos1]))/2
            if 0<auc:
                peaks.append(Peak(peak_loc.rt,peak_loc.sc,len(rd),peak_loc.coef,auc,rt_all[pos0],rt_all[pos1-1]))


    else:
        for rd in ridgelines:
            peak_loc=max(rd,key=operator.attrgetter('coef'))

            pos=bisect_left(rt_all,peak_loc.rt)
            auc=eic_dict[rt_all[pos-1]]+eic_dict[rt_all[pos]]
            peaks.append(Peak(peak_loc.rt,peak_loc.sc,len(rd),peak_loc.coef,auc,None,None))






    peaks.sort(key=operator.attrgetter('coef'),reverse=True)
    peaks_=peaks[:]
    for nn,peak1 in enumerate(peaks_):
        for p in peaks_[nn+1:]:
            if p in peaks and abs(peak1.rt-p.rt)<.1:
                peaks.remove(p)
    if not peaks:
        sc=10
        midpt=(rt_all[0]+rt_all[-1])/2
        return [Peak(midpt,0,0,0,1,midpt-sc,midpt+sc)]


    return[x for x in peaks if x.auc>peaks[0].auc/10]
    return peaks

