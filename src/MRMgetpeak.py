import collections
import glob
import operator
from functools import partial
import sys
import statistics
import math
import os
import concurrent.futures
import re
from bisect import bisect_left
from bisect import bisect_right
from multiprocessing import freeze_support


import commonfn
import MRMcwt

param_set = {
    "mzML_files",
    "transition_list",
    "ISTD_trace_all",
    "max_RTshift",
    "normalise",
}
param_dict = commonfn.read_param(param_set)

normalise_bool = int(param_dict.get("normalise", 1))

mzML_files = []
for x in param_dict["mzML_files"].split("\t"):
    mzML_files.extend(glob.glob(x))
mzML_files.sort()

transition = param_dict["transition_list"]
ISTD_trace_all = int(param_dict["ISTD_trace_all"])

t_list = commonfn.read_trans_stat(transition)
istd_dict = commonfn.get_istd_dict(t_list)

miscdir = "trans"
trans2filename = {
    trans: commonfn.format_filename(
        "{}Q1_{:.1f}Q3_{:.1f}.txt".format(trans.name, trans.q1, trans.q3)
    )
    for trans in t_list
}  # +t_list2}
trans1filename = {
    trans: os.path.join(miscdir, "trans_" + t) for trans, t in trans2filename.items()
}
trans2filename = {
    trans: os.path.join(miscdir, "trans2_" + t) for trans, t in trans2filename.items()
}

t_list = [
    x for x in t_list if x.istd in istd_dict and os.path.isfile(trans1filename[x])
]


wave_scales = [2, 2.5] + list(range(3, 9))
wave_sqrt = [math.sqrt(x) for x in wave_scales]
findridge = partial(MRMcwt.findridge0, wave_scales, wave_sqrt, True)


def get_EIC_feats(feat_in):
    next(feat_in)  # name
    RTline = feat_in.readline()  # RT line
    while RTline:
        count_feats = 0
        rt_l = [float(x) for x in RTline.split("\t")]
        I_l = [float(x) for x in feat_in.readline().split("\t")]
        feat0 = feat_in.readline().strip()
        feat_dat = []
        while feat0:
            count_feats += 1
            if count_feats == 1:
                feat_dat.append([float(x) for x in feat0.split("\t")])
            feat0 = feat_in.readline().strip()
        yield rt_l, I_l, feat_dat
        RTline = feat_in.readline()


max_RTshift = float(param_dict.get("max_RTshift", 5))


def my_median(ll):
    lll = sorted(ll)
    return lll[int(len(lll) / 2)]


def write_to_trans(trans):
    istd0 = istd_dict[trans.istd]
    if istd0.istd != istd0.name:
        istd0 = istd_dict[istd0.istd]
    if not os.path.isfile(trans1filename[istd0]):
        return

    with open(trans1filename[trans]) as feat_in, open(trans2filename[trans], "w") as fo:

        print(trans.name)

        trans_feats = []
        for eicfeat in get_EIC_feats(feat_in):
            rt_, I_, *_ = eicfeat
            firstdecile = sorted(I_)[int(len(I_) / 10)]
            I_ = [max(0, i - firstdecile) for i in I_]
            rerunw = 2 * (rt_[-1] - rt_[0]) / (len(rt_) - 1)
            feats = findridge(rerunw, list(zip(rt_, I_)))
            trans_feats.append(feats)
        feat_in.seek(0)

        if trans.name in istd_dict:
            med_npeaks = 1
        else:
            med_npeaks = sorted(len(x) for x in trans_feats)[int(len(trans_feats) / 2)]

        median_rt = []
        for npeaks in range(1, med_npeaks + 1):
            feats_rt = []
            for feats in (x for x in trans_feats if npeaks <= len(x)):
                if npeaks < len(x):
                    feats = sorted(feats, key=operator.attrgetter("auc"))[-npeaks:]
                feats_rt.append(sorted(x.rt for x in feats))
            cur_median_rt = [my_median(rt[i] for rt in feats_rt) for i in range(npeaks)]
            if all(any(abs(x - y) < 3 for y in cur_median_rt) for x in median_rt):
                median_rt = cur_median_rt
            else:
                break

        if not median_rt:
            median_rt = [(rt_[0] + rt_[-1]) / 2]

        istd_in = open(trans1filename[istd0])
        fo.write("{}\n".format(trans.name))
        for ii, (eicfeat, istd_eicfeat) in enumerate(
            zip(get_EIC_feats(feat_in), get_EIC_feats(istd_in))
        ):
            rt_, I_, *_ = eicfeat
            firstdecile = sorted(I_)[int(len(I_) / 10)]
            I_ = [max(0, i - firstdecile) for i in I_]

            maxf = istd_eicfeat[2][0]

            feats_rt = list(x.rt for x in trans_feats[ii])
            new_rt = []
            added = set()
            for feat_rt in median_rt:
                remaining = [x for x in feats_rt if x not in added]
                if remaining:
                    feat0 = min(remaining, key=lambda x: abs(x - feat_rt))
                    if abs(feat0 - feat_rt) < max_RTshift:  # or len(median_rt)==1:
                        new_rt.append((feat0, 1))
                        added.add(feat0)
                    else:
                        new_rt.append((feat_rt, 3))
                else:
                    new_rt.append((feat_rt, 3))
            new_rt.sort()

            fo.write("{}\n".format("\t".join(str(x) for x in rt_)))
            fo.write("{}\n".format("\t".join(str(x) for x in eicfeat[1])))

            firstdecile_ = sorted(istd_eicfeat[1])[int(len(istd_eicfeat[1]) / 10)]
            RT_I = [
                (x, max(0, y - firstdecile_))
                for x, y in zip(*istd_eicfeat[:2])
                if maxf[2] - 0.001 < x < maxf[3] + 0.001
            ]
            if 2 * (maxf[0] - maxf[2]) < maxf[3] - maxf[0]:
                RT_I = [
                    (2 * maxf[0] - x, y)
                    for x, y in RT_I[::-1]
                    if x > 1 + 2 * maxf[0] - maxf[2]
                ] + RT_I
            elif 2 * (maxf[3] - maxf[0]) < maxf[0] - maxf[2]:
                RT_I = RT_I + [
                    (2 * maxf[0] - x, y)
                    for x, y in RT_I[::-1]
                    if x < -1 + 2 * maxf[0] - maxf[3]
                ]
            topN = sorted(RT_I, key=lambda x: abs(x[0] - maxf[0]))[:4]
            apex_istd = max(x for _, x in topN + [(None, 1)])  # min 1

            if ISTD_trace_all:
                for feat_rt, cc in new_rt:
                    topN = sorted(zip(rt_, I_), key=lambda x: abs(x[0] - feat_rt))[:2]
                    apex_int = max(x for _, x in topN)
                    RT_I_ = [(rt, x / apex_istd * apex_int) for rt, x in RT_I]
                    auc = (
                        sum(
                            (rt_i0[1] + rt_i1[1]) * (rt_i1[0] - rt_i0[0])
                            for rt_i0, rt_i1 in zip(RT_I_, RT_I_[1:])
                        )
                        / 2
                    )
                    fo.write("{:.3f}\t{}\t{}\t".format(feat_rt, auc, cc))
                    fo.write("\t".join(str(x - maxf[0] + feat_rt) for x, _ in RT_I))
                    fo.write("\t")
                    fo.write(
                        "\t".join(
                            str(x / apex_istd * apex_int + firstdecile) for _, x in RT_I
                        )
                    )
                    fo.write("\n")
                fo.write("\n")
            else:
                for jj, (feat_rt, cc) in enumerate(new_rt):
                    if jj == 0:
                        pos0 = bisect_left(rt_, feat_rt - 2.5 * maxf[4])
                    else:
                        pos0 = bisect_left(rt_, (new_rt[jj - 1][0] + feat_rt) / 2)
                    if jj == len(new_rt) - 1:
                        pos1 = bisect_left(rt_, feat_rt + 2.5 * maxf[4])
                    else:
                        pos1 = bisect_left(rt_, (new_rt[jj + 1][0] + feat_rt) / 2)
                    RT_I_ = list(zip(rt_, I_))[pos0:pos1]
                    auc = (
                        sum(
                            (rt_i0[1] + rt_i1[1]) * (rt_i1[0] - rt_i0[0])
                            for rt_i0, rt_i1 in zip(RT_I_, RT_I_[1:])
                        )
                        / 2
                    )
                    fo.write("{:.3f}\t{}\t{}\t".format(feat_rt, auc, cc))

                    fo.write("{}\t".format("\t".join(str(x) for x in rt_[pos0:pos1])))
                    fo.write(
                        "{}\n".format("\t".join(str(x) for x in eicfeat[1][pos0:pos1]))
                    )

                fo.write("\n")


if __name__ == "__main__":
    print(len(mzML_files), "mzML files")

    ppfile = "peak_picking.txt"
    if os.path.isfile(ppfile):
        sorted_name = [x.name for x in t_list]
        t_list_sub = dict()
        print(ppfile + " exists, re-integrating")
        with open(ppfile) as usersRT:
            next(usersRT)
            for line in usersRT:
                lsp = line.rstrip("\n").split("\t")
                rt_l = re.findall("\d+\.*\d*-\d+\.*\d*", lsp[8])
                if rt_l:
                    pos0 = bisect_left(sorted_name, lsp[1])
                    pos1 = bisect_right(sorted_name, lsp[1])
                    for trans in t_list[pos0:pos1]:
                        if (trans.q1, trans.q3) == tuple(float(x) for x in lsp[2:4]):
                            t_list_sub[trans] = [float(x) for x in rt_l[0].split("-")]
                            break
        for trans, (beg0, end0) in t_list_sub.items():
            with open(trans1filename[trans]) as feat_in, open(
                trans2filename[trans], "w"
            ) as fo:
                print(trans.name)
                fo.write("{}\n".format(trans.name))
                for eicfeat in get_EIC_feats(feat_in):
                    rt_, I_, *_ = eicfeat
                    firstdecile = sorted(I_)[int(len(I_) / 10)]
                    I_ = [max(0, i - firstdecile) for i in I_]
                    pos0 = bisect_left(rt_, beg0)
                    pos1 = bisect_left(rt_, end0)
                    RT_I_ = list(zip(rt_, I_))[pos0:pos1]
                    auc = (
                        sum(
                            (rt_i0[1] + rt_i1[1]) * (rt_i1[0] - rt_i0[0])
                            for rt_i0, rt_i1 in zip(RT_I_, RT_I_[1:])
                        )
                        / 2
                    )
                    fo.write("{}\n".format("\t".join(str(x) for x in rt_)))
                    fo.write("{}\n".format("\t".join(str(x) for x in eicfeat[1])))
                    fo.write("{:.3f}\t{}\t1\t".format((beg0 + end0) / 2, auc))
                    fo.write("{}\t".format("\t".join(str(x) for x in rt_[pos0:pos1])))
                    fo.write(
                        "{}\n\n".format(
                            "\t".join(str(x) for x in eicfeat[1][pos0:pos1])
                        )
                    )

    else:
        freeze_support()
        with concurrent.futures.ProcessPoolExecutor(max_workers=19) as executor:
            list(executor.map(write_to_trans, t_list))

    mzML_files = [mf.rstrip() for mf in open("trans_mzML_list.txt")]

    def qq_rt(secondr):
        wstr0 = "trans2_" if secondr else "trans_"
        QQ_rt = collections.defaultdict(list)
        for trans in t_list:
            trans_file = os.path.join(
                miscdir,
                commonfn.format_filename(
                    wstr0
                    + "{}Q1_{:.1f}Q3_{:.1f}.txt".format(trans.name, trans.q1, trans.q3)
                ),
            )
            if os.path.isfile(trans_file):
                fo = open(trans_file)
                next(fo)
                while fo.readline():
                    next(fo)
                    feat0 = fo.readline().strip()
                    feat_l = []
                    while feat0:
                        feat_l.append([float(x) for x in feat0.split("\t")[:2]])
                        feat0 = fo.readline().strip()
                    QQ_rt[trans].append(feat_l)
        return dict(QQ_rt)

    QQ_rt2 = qq_rt(1)

    with open("quant_auc_rt.txt", "w") as quant_:
        quant_.write(
            "name\tprecursor\tproduct\tRT\t{}\t{}\t{}\n".format(
                "\t".join(mzML_files),
                "\t".join("norm_" + x for x in mzML_files),
                "\t".join("npeaks_" + x for x in mzML_files),
            )
        )
        for trans, auc_rt_l in QQ_rt2.items():

            istd0 = istd_dict[trans.istd]
            for jj in range(len(auc_rt_l[0])):
                medianRT = statistics.median(
                    auc_rt_l[ii][jj][0] for ii in range(len(mzML_files))
                )
                quant_.write(
                    "{} RT{:.0f}\t{:.1f}\t{:.1f}\t{}".format(
                        trans.name, medianRT, trans.q1, trans.q3, trans.rt
                    )
                )
                for ii in range(len(mzML_files)):
                    quant_.write("\t{}".format(auc_rt_l[ii][jj][1]))

                if normalise_bool:
                    if trans.name in istd_dict:
                        quant_.write("\t" * len(mzML_files))
                    else:
                        auc_istd = QQ_rt2[istd0]
                        for ii in range(len(mzML_files)):
                            auc_istd_ = auc_istd[ii][0][1]
                            if auc_istd_ == 0:
                                quant_.write("\t{}".format(0))
                            else:
                                quant_.write(
                                    "\t{}".format(auc_rt_l[ii][jj][1] / auc_istd_)
                                )
                else:
                    for ii in range(len(mzML_files)):
                        quant_.write("\t{}".format(auc_rt_l[ii][jj][1]))

                for ii in range(len(mzML_files)):
                    quant_.write("\t{}".format(len(auc_rt_l[ii])))

                quant_.write("\n")
