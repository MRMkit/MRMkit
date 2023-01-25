import collections
import operator
from bisect import bisect_left
import re
import os

import commonfn


t_list = dict()
with open("peak_picking.txt") as usersRT:
    next(usersRT)
    for line in usersRT:
        lsp = line.rstrip("\n").split("\t")
        rt_l = re.findall("\d+\.*\d*-\d+\.*\d*", lsp[9])
        if rt_l:
            t_list[
                "{}Q1_{:.1f}Q3_{:.1f}".format(lsp[1], float(lsp[2]), float(lsp[3]))
            ] = [sum(float(x) for x in rt_l[0].split("-")) / 2]
        else:
            rt_l = re.findall("\d+\.*\d*", lsp[9])
            if rt_l:
                t_list[
                    "{}Q1_{:.1f}Q3_{:.1f}".format(lsp[1], float(lsp[2]), float(lsp[3]))
                ] = [float(x) for x in rt_l]


QCdict = dict()
with open("BQC_CoV.txt") as BQCtab, open("TQC_CoV.txt") as TQCtab, open(
    "D_ratio.txt"
) as D_ratio:
    BQCtab.readline()
    TQCtab.readline()
    D_ratio.readline()
    for bline, tline, rline in zip(BQCtab, TQCtab, D_ratio):
        blsp = bline.rstrip("\n").split("\t")[:6]
        tlsp = tline.rstrip("\n").split("\t")[:6]
        rlsp = rline.rstrip("\n").split("\t")[:6]
        QCdict[
            (blsp[0].split(" RT")[0],)
            + tuple(blsp[1:4])
            + (int(blsp[0].split(" RT")[1]),)
        ] = (blsp[5], tlsp[5], rlsp[5])


with open("batch_adjusted.txt") as batchadj, open("quant_table.txt", "w") as batchadj_:
    cpdname = batchadj.readline().rstrip("\n").split("\t")[3:]
    q1line = batchadj.readline().rstrip().split("\t")[3:]
    q3line = batchadj.readline().rstrip().split("\t")[3:]
    RTline = batchadj.readline().rstrip().split("\t")[3:]
    linesp = [line.rstrip("\n").split("\t")[3:] for line in batchadj]
    batchadj.seek(0)
    first3 = ("\t".join(line.split("\t", 3)[:3]) for line in batchadj)

    cpdRT = collections.defaultdict(list)
    col_dict = dict()
    for ii, (cpd_RT, q1, q3, rt) in enumerate(zip(cpdname, q1line, q3line, RTline)):
        col_dict[(cpd_RT, q1, q3, rt)] = [lsp[ii] for lsp in linesp]
        cpd, RT = cpd_RT.split(" RT", 1)
        cpdRT[(cpd, q1, q3, rt)].append(int(RT))
    selcpd = dict()
    for (cpd, q1, q3, rt), rt_l in cpdRT.items():
        user_rt = t_list.get("{}Q1_{:.1f}Q3_{:.1f}".format(cpd, float(q1), float(q3)))
        if user_rt is not None:
            rt_ = [min(rt_l, key=lambda x: abs(x - y)) for y in user_rt]
            selcpd[(cpd, q1, q3, rt)] = rt_
    sumcol_dict = dict()
    for (cpd, q1, q3, rt), rt_l in selcpd.items():
        colv = [float(x) for x in col_dict[(cpd + " RT" + str(rt_l[0]), q1, q3, rt)]]
        for rt_ in rt_l[1:]:
            colv = [
                x + float(y)
                for x, y in zip(colv, col_dict[(cpd + " RT" + str(rt_), q1, q3, rt)])
            ]
        sumcol_dict[(cpd, q1, q3, rt)] = colv

    batchadj_.write(next(first3) + "\t" + "\t".join(x for x, _, _, _ in selcpd) + "\n")
    batchadj_.write(next(first3) + "\t" + "\t".join(x for _, x, _, _ in selcpd) + "\n")
    batchadj_.write(next(first3) + "\t" + "\t".join(x for _, _, x, _ in selcpd) + "\n")
    next(first3)
    batchadj_.write(
        "detected_RT\t\t"
        + "\t"
        + "\t".join(", ".join(str(xx) for xx in x) for x in selcpd.values())
        + "\n"
    )
    batchadj_.write(
        "%CoV(BQC)\t\t"
        + "\t"
        + "\t".join(
            ", ".join(QCdict[k + (xx,)][0] for xx in x) for k, x in selcpd.items()
        )
        + "\n"
    )
    batchadj_.write(
        "%CoV(TQC)\t\t"
        + "\t"
        + "\t".join(
            ", ".join(QCdict[k + (xx,)][1] for xx in x) for k, x in selcpd.items()
        )
        + "\n"
    )
    batchadj_.write(
        "D-ratio\t\t"
        + "\t"
        + "\t".join(
            ", ".join(QCdict[k + (xx,)][2] for xx in x) for k, x in selcpd.items()
        )
        + "\n"
    )
    for i in range(len(colv)):
        batchadj_.write(next(first3))
        for (cpd, q1, q3, rt), colv in sumcol_dict.items():
            batchadj_.write("\t{}".format(colv[i]))
        batchadj_.write("\n")
