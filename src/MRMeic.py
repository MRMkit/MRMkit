import xml.etree.ElementTree as ET
import base64
import struct
import zlib

import commonfn

Eic = commonfn.Eic


def bin2float(node):
    d = base64.b64decode(node.findtext("{http://psi.hupo.org/ms/mzml}binary"))
    if node.find("*/[@accession='MS:1000574']") is not None:
        d = zlib.decompress(d)
    fmt = (
        "<" + str(int(len(d) / 4)) + "f"
        if node.find("*/[@accession='MS:1000523']") is None
        else "<" + str(int(len(d) / 8)) + "d"
    )
    return struct.unpack(fmt, d)


def store_scan(element):
    time_arr = bin2float(element.find(".//*[@accession='MS:1000595'].."))
    if element.find(".//*[@accession='MS:1000595']").attrib["unitName"] == "minute":
        time_arr = [x * 60 for x in time_arr]
    inten = bin2float(element.find(".//*[@accession='MS:1000515'].."))
    return time_arr, inten
    return [m for m, i in zip(time_arr, inten) if i > 0], [i for i in inten if i > 0]


def print_eic(mzML_file):

    eic_dict0 = dict()
    eickey = 0

    for _, element in ET.iterparse(mzML_file):
        if element.tag == "{http://psi.hupo.org/ms/mzml}chromatogram":
            precursor = element.find(
                "{http://psi.hupo.org/ms/mzml}precursor//*[@accession='MS:1000827']"
            )
            product = element.find(
                "{http://psi.hupo.org/ms/mzml}product//*[@accession='MS:1000827']"
            )
            Q1, Q3 = 0, 0
            if precursor is not None and product is not None:
                Q1, Q3 = float(precursor.attrib["value"]), float(
                    product.attrib["value"]
                )
            rt_, I_ = store_scan(element)
            eic_dict0[Eic(Q1, Q3, eickey, min(rt_), max(rt_))] = (rt_, I_)
            eickey += 1
    return {k: eic_dict0[k] for k in sorted(eic_dict0)}
