import collections
import operator
from bisect import bisect_left
import re


t_list=dict()
with open('peak_picking.txt') as usersRT:
    next(usersRT)
    for line in usersRT:
        lsp=line.rstrip('\n').split('\t')
        rt_l=re.findall("\d+\.*\d*",lsp[8])
        if rt_l:
            t_list['{}\t{:.1f}\t{:.1f}'.format(lsp[1],float(lsp[2]),float(lsp[3]))]=[float(x)for x in rt_l]



with open('batch_adjusted.txt') as batchadj, open('quant_table.txt','w') as batchadj_:
    cpdname=batchadj.readline().rstrip('\n').split('\t')[3:]
    q1line=batchadj.readline().rstrip().split('\t')[3:]
    q3line=batchadj.readline().rstrip().split('\t')[3:]
    RTline=batchadj.readline().rstrip().split('\t')[3:]
    linesp=[line.rstrip('\n').split('\t')[3:]for line in batchadj]
    batchadj.seek(0)
    first3=('\t'.join(line.split('\t',3)[:3]) for line in batchadj)

    cpdRT=collections.defaultdict(list)
    col_dict=dict()
    for ii,(cpd_RT,q1,q3,rt) in enumerate(zip(cpdname,q1line,q3line,RTline)):
        col_dict[(cpd_RT,q1,q3,rt)]=[lsp[ii] for lsp in linesp]
        cpd,RT=cpd_RT.split(' RT',1)
        cpdRT[(cpd,q1,q3,rt)].append(int(RT))
    selcpd=dict()
    for (cpd,q1,q3,rt),rt_l in cpdRT.items():
        user_rt=t_list.get('{}\t{:.1f}\t{:.1f}'.format(cpd,float(q1),float(q3)))
        if user_rt is not None:
            rt_=[min(rt_l,key=lambda x:abs(x-y))for y in user_rt]
            selcpd[(cpd,q1,q3,rt)]=rt_
    sumcol_dict=dict()
    for (cpd,q1,q3,rt),rt_l in selcpd.items():
        colv=[float(x) for x in col_dict[(cpd+' RT'+str(rt_l[0]),q1,q3,rt)]]
        for rt_ in rt_l[1:]:
            colv=[x+float(y) for x,y in zip(colv,col_dict[(cpd+' RT'+str(rt_),q1,q3,rt)])]
        sumcol_dict[(cpd,q1,q3,rt)]=colv

    batchadj_.write(next(first3)+'\t'+'\t'.join(x for x,_,_,_ in selcpd)+'\n')
    batchadj_.write(next(first3)+'\t'+'\t'.join(x for _,x,_,_ in selcpd)+'\n')
    batchadj_.write(next(first3)+'\t'+'\t'.join(x for _,_,x,_ in selcpd)+'\n')
    batchadj_.write(next(first3)+'\t'+'\t'.join(x for _,_,_,x in selcpd)+'\n')
    for i in range(len(colv)):
        batchadj_.write(next(first3))
        for (cpd,q1,q3,rt),colv in sumcol_dict.items():
            batchadj_.write('\t{}'.format(colv[i]))
        batchadj_.write('\n')
