from bisect import bisect_left
import operator
import os
import collections
import sys

Peak = collections.namedtuple("Peak", ("rt sc rlen coef auc begin end"))
Eic = collections.namedtuple("Eic", ("q1 q3 eickey minRT maxRT"))
QQ = collections.namedtuple("QQ", ("lclass q1 q3 name rt istd"))


def read_param(param_set):
    param_dict = dict()
    with open("param.txt") as input_param:
        key_ = ""
        value_ = ""
        for line in (l.strip() for l in input_param if l[0] != "#"):
            if not line and key_ and value_:
                param_dict[key_] = value_[:-1]
                key_ = ""
                value_ = ""
            elif line in param_set:
                key_ = line
            elif key_:
                value_ += line + "\n"
        param_dict[key_] = value_[:-1]
    return param_dict


def read_trans_stat(transition):
    t_list = set()
    with open(transition) as trans:
        trans.readline()
        for line in (l for l in trans if l.strip() and l[0] != "#"):
            lsp = [x.strip() for x in line.split("\t")]
            if not lsp[2].strip():
                lsp[2] = lsp[1]
            t_list.add(
                QQ(lsp[0], float(lsp[3]), float(lsp[4]), lsp[1], float(lsp[5]), lsp[2])
            )
    return sorted(t_list, key=operator.attrgetter("name"))


def get_istd_dict(t_list):
    istd_dict = dict()
    t_list.sort(key=operator.attrgetter("name"))
    sorted_name = [x.name for x in t_list]
    for trans in t_list:
        if trans.istd not in istd_dict:
            pos0 = bisect_left(sorted_name, trans.istd)
            if t_list[pos0].name == trans.istd:
                istd_dict[trans.istd] = t_list[pos0]
            else:
                print("ISTD not found for {}".format(trans.name))
    return istd_dict


def find_trans(Q1Q3, trans):
    pos0 = bisect_left(Q1Q3, (trans.q1 - 0.09,))
    pos1 = bisect_left(Q1Q3, (trans.q1 + 0.09,))
    filter0 = [
        x
        for x in Q1Q3[pos0:pos1]
        if abs(x.q3 - trans.q3) < 0.09 and x.minRT < trans.rt < x.maxRT
    ]
    if filter0:
        return min(filter0, key=lambda x: abs((x.minRT + x.maxRT) / 2 - trans.rt))
    return None


import string


def format_filename(s):
    valid_chars = "-_.(){}{}".format(string.ascii_letters, string.digits)
    filename = "".join((c if c in valid_chars else "_") for c in s)
    return filename
