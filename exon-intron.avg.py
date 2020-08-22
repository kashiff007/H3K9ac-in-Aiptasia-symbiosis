:q!#!/usr/bin/python

import sys
import numpy as np
np.set_printoptions(threshold=np.inf)
#from itertools import count, tee, izip, islice
import GFFutils
import pyBigWig
import matplotlib.pyplot as plt

G = GFFutils.GFFDB("/home/user/dm1.db")
C_bw = pyBigWig.open("/home/user/Symb_treatedVscontrol_50bin.bw")
Gene_final = []
def separate_exon_intron(EI_list, EI1, EI2, EI3):
    for idx, ele in enumerate(EI_list):
        ele = len(ele)
        if int(idx) == 0:
            EI1.append(ele)
        elif int(idx) == 1:
            EI2.append(ele)
        elif int(idx) == 2:
            EI3.append(ele)

exon_len = []
upstream_3000, exon1, intron1, exon2, intron2, exon3,intron3, after_1000, before_1000,intron_3,exon_3, intron_2, exon_2, intron_1, exon_1, downstream_3000 = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]


for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))

    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    separate_exon_intron(first_3_exons, exon1, exon2, exon3)
    separate_exon_intron(last_3_exons, exon_3, exon_2, exon_1)
    separate_exon_intron(first_3_introns, intron1, intron2, intron3)
    separate_exon_intron(last_3_introns, intron_3, intron_2, intron_1)

exon1_mean,exon2_mean,exon3_mean,exon_1_mean,exon_2_mean,exon_3_mean,intron1_mean,intron2_mean,intron3_mean,intron_1_mean,intron_2_mean,intron_3_mean = round(np.mean(exon1),2), round(np.mean(exon2),2), round(np.mean(exon3),2), round(np.mean(exon_1),2), round(np.mean(exon_2),2), round(np.mean(exon_3),2), round(np.mean(intron1),2), round(np.mean(intron2),2), round(np.mean(intron3),2), round(np.mean(intron_1),2), round(np.mean(intron_2),2), round(np.mean(intron_3),2)
    
print("exon1: ", exon1_mean, "frequency: ", len(exon1))
print("exon2: ", exon2_mean, "frequency: ", len(exon2))
print("exon3: ", exon3_mean, "frequency: ", len(exon3))
print("exon_3: ", exon_3_mean, "frequency: ", len(exon_3))
print("exon_2: ", exon_2_mean, "frequency: ", len(exon_2))
print("exon_1: ", exon_1_mean, "frequency: ", len(exon_1))

print("intron1: ", intron1_mean, "frequency: ", len(intron1))
print("intron2: ", intron2_mean, "frequency: ", len(intron2))
print("intron3: ", intron3_mean, "frequency: ", len(intron3))
print("intron_3: ", intron_3_mean, "frequency: ", len(intron_3))
print("intron_2: ", intron_2_mean, "frequency: ", len(intron_2))
print("intron_1: ", intron_1_mean, "frequency: ", len(intron_1))

def divide_chunks(l, n): 
      
    # looping till length l 
    if len(l) > 9:
        for i in range(0, len(l), n):
            yield (l[i:i + n])

'''  
def normalized_region_score(region, num, EI):
    
    final1 =[]
    for mRNA in G.features_of_type('mRNA'):
        exons = list(G.children(mRNA, featuretype='exon'))
        introns = list(G.interfeatures(exons))
        if mRNA.strand == "-":
            first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
        else:
            first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

        for idx, ele in enumerate(region):
            if int(idx) == num:
                if ele.chrom in C_bw.chroms():
                    vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                    exon_len = len(vals) - len(vals) % 10
                    vals = vals[:exon_len]
                    bin_num = (exon_len/10)
                    x = (list(divide_chunks(vals, bin_num)))
                    if x == []:
                        continue
                    a = np.array(x)
                    y = np.mean(a, axis=1)
                    final1.append(y.tolist())
    final1_mean = np.array(final1)
    print ("exon 1: ", np.nanmean(final1_mean, axis=0))
    EI = np.nanmean(final1_mean, axis=0)
    print EI
normalized_region_score(first_3_exons, 0, exon1)
'''

#upstream 3000bp
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_exons):
        if int(idx) == 0:
            if ele.chrom in C_bw.chroms():
                if mRNA.strand == "-":
                    v = ele.stop + 3000
                    if (int(v) > int(C_bw.chroms(ele.chrom))) and (ele.stop - int(C_bw.chroms(ele.chrom))):
                        x = [0]*((v) - int(C_bw.chroms(ele.chrom)))
                        vals = C_bw.values(ele.chrom, ele.stop, int(C_bw.chroms(ele.chrom)))
                        a = x + vals
                        x = (list(divide_chunks(a, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    elif int(C_bw.chroms(ele.chrom)) > v:
                        vals = C_bw.values(ele.chrom, ele.stop, v)
                        x = (list(divide_chunks(vals, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                else:
                    v = ele.start - 3000
                    if v < 0:
                        x = [0] * int(3000-(ele.start))
                        vals = C_bw.values(ele.chrom, 0, int(ele.start))
                        a = x + vals
                        x = (list(divide_chunks(a, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    else:
                        vals = C_bw.values(ele.chrom, v, ele.start)
                        x = (list(divide_chunks(vals, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
            else:
                x = [0]*50
                a = np.array(x)
                final1.append(a.tolist())


final1_mean = np.array(final1)
print ("upstream 3000bp: ", np.nanmean(final1_mean, axis=0))
upstream_3000bp = (np.nanmean(final1_mean, axis=0))

#exon1     
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_exons):
        if int(idx) == 0:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = (exon_len/10)
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("exon 1: ", np.nanmean(final1_mean, axis=0))
exon1 = np.nanmean(final1_mean, axis=0)

#intron1
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_introns):
        if int(idx) == 0:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("intron 1: ", np.nanmean(final1_mean, axis=0))
intron1 = np.nanmean(final1_mean, axis=0)

#exon2
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_exons):
        if int(idx) == 1:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("exon 2: ", np.nanmean(final1_mean, axis=0))
exon2 = np.nanmean(final1_mean, axis=0)

#intron2
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_introns):
        if int(idx) == 1:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())

final1_mean = np.array(final1)
print ("intron 2: ", np.nanmean(final1_mean, axis=0))
intron2 = np.nanmean(final1_mean, axis=0)

#exon3
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_exons):
        if int(idx) == 2:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("exon 3: ", np.nanmean(final1_mean, axis=0))
exon3 = np.nanmean(final1_mean, axis=0)

#after_1000
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(first_3_exons):
        if int(idx) == 2:
            if ele.chrom in C_bw.chroms():
                if mRNA.strand == "-":
                    v = ele.start - 1000
                    if int(v) < 0:
                        x = [0]*((1000) - ele.start)
                        vals = C_bw.values(ele.chrom, 0, ele.start)
                        a = x + vals
                        x = (list(divide_chunks(a, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    else:
                        vals = C_bw.values(ele.chrom, v, ele.start)
                        x = (list(divide_chunks(vals, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                else:
                    v = ele.stop + 1000
                    if (v > int(C_bw.chroms(ele.chrom)) and (ele.stop - int(C_bw.chroms(ele.chrom)))): 
                        x = [0] * int(v-int(C_bw.chroms(ele.chrom)))
                        vals = C_bw.values(ele.chrom, ele.stop, int(C_bw.chroms(ele.chrom)))
                        a = x + vals
                        x = (list(divide_chunks(a, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    elif int(C_bw.chroms(ele.chrom)) > v:
                        vals = C_bw.values(ele.chrom, ele.stop, v)
                        x = (list(divide_chunks(vals, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())


final1_mean = np.array(final1)
print ("after 1000bp: ", np.nanmean(final1_mean, axis=0))
after_1000 = np.nanmean(final1_mean, axis=0)

#before_1000
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_exons):
        if int(idx) == 0:
            if ele.chrom in C_bw.chroms():
                if mRNA.strand == "-":
                    v = ele.stop + 1000
                    if (int(v) > int(C_bw.chroms(ele.chrom)) and (ele.stop - int(C_bw.chroms(ele.chrom)))):
                        x = [0]*(v - int(C_bw.chroms(ele.chrom)))
                        vals = C_bw.values(ele.chrom, ele.stop, int(C_bw.chroms(ele.chrom)))
                        a = x + vals
                        x = (list(divide_chunks(a, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    elif int(C_bw.chroms(ele.chrom)) > v:
                        vals = C_bw.values(ele.chrom, ele.stop, v)
                        x = (list(divide_chunks(vals, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                else:
                    v = ele.start - 1000
                    if v < 0: 
                        x = [0] * int(1000-ele.start)
                        vals = C_bw.values(ele.chrom, 0, ele.start)
                        a = x + vals
                        x = (list(divide_chunks(a, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    else:
                        vals = C_bw.values(ele.chrom, v, ele.start)
                        x = (list(divide_chunks(vals, 100)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())


final1_mean = np.array(final1)
print ("Before exon -3 1000bp", np.nanmean(final1_mean, axis=0))
before_1000 = np.nanmean(final1_mean, axis=0)


#exon-3     
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_exons):
        if int(idx) == 0:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("exon -3: ", np.nanmean(final1_mean, axis=0))
exon_3 = np.nanmean(final1_mean, axis=0)

#intron-2
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_introns):
        if int(idx) == 1:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("intron -2: ", np.nanmean(final1_mean, axis=0))
intron_2 = np.nanmean(final1_mean, axis=0)

#exon-2 
final1 = []
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_exons):
        if int(idx) == 1:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())

final1_mean = np.array(final1)
print ("exon -2: ", np.nanmean(final1_mean, axis=0))
exon_2 = np.nanmean(final1_mean, axis=0)

#intron-1
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_introns):
        if int(idx) == 2:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("intron -1: ", np.nanmean(final1_mean, axis=0))
intron_1 = np.nanmean(final1_mean, axis=0)

#exon-1     
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_exons):
        if int(idx) == 2:
            if ele.chrom in C_bw.chroms():
                vals = C_bw.values(ele.chrom, ele.start, ele.stop)
                exon_len = len(vals) - len(vals) % 10
                vals = vals[:exon_len]
                bin_num = exon_len/10
                x = (list(divide_chunks(vals, bin_num)))
                if x == []:
                    continue
                a = np.array(x)
                y = np.mean(a, axis=1)
                final1.append(y.tolist())
            else:
                x = [0]*10
                a = np.array(x)
                final1.append(a.tolist())
final1_mean = np.array(final1)
print ("exon -1: ", np.nanmean(final1_mean, axis=0))
exon_1 = np.nanmean(final1_mean, axis=0)


#downstream 3000bp
final1 =[]
for mRNA in G.features_of_type('mRNA'):
    exons = list(G.children(mRNA, featuretype='exon'))
    introns = list(G.interfeatures(exons))
    if mRNA.strand == "-":
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[-3:], exons[:3], introns[-3:], introns[:3]
    else:
        first_3_exons,last_3_exons,first_3_introns,last_3_introns = exons[:3], exons[-3:], introns[:3], introns[-3:]

    for idx, ele in enumerate(last_3_exons):
        if int(idx) == 0:
            if ele.chrom in C_bw.chroms():
                if mRNA.strand == "+":
                    v = ele.stop + 3000
                    if (int(v) > int(C_bw.chroms(ele.chrom))) and (ele.stop - int(C_bw.chroms(ele.chrom))):
                        x = [0]*((v) - int(C_bw.chroms(ele.chrom)))
                        vals = C_bw.values(ele.chrom, ele.stop, int(C_bw.chroms(ele.chrom)))
                        a = x + vals
                        x = (list(divide_chunks(a, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    elif int(C_bw.chroms(ele.chrom)) > v:
                        vals = C_bw.values(ele.chrom, ele.stop, v)
                        x = (list(divide_chunks(vals, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                else:
                    v = ele.start - 3000
                    if v < 0:
                        x = [0] * int(3000-(ele.start))
                        vals = C_bw.values(ele.chrom, 0, int(ele.start))
                        a = x + vals
                        x = (list(divide_chunks(a, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
                    else:
                        vals = C_bw.values(ele.chrom, v, ele.start)
                        x = (list(divide_chunks(vals, 60)))
                        if x == []:
                            continue
                        y = np.mean(np.array(x), axis=1)
                        final1.append(y.tolist())
            else:
                x = [0]*50
                a = np.array(x)
                final1.append(a.tolist())


final1_mean = np.array(final1)
print ("downstream 3000bp: ", np.nanmean(final1_mean, axis=0))
downstream_3000 = np.nanmean(final1_mean, axis=0)

#Gene_len1 =  np.concatenate((upstream_3000 , exon1 , intron1 , exon2 , intron2 , exon3 , after_1000 , before_1000 , exon_3 , intron_2 , exon_2 , intron_1 , exon_1 , downstream_3000), axis=None)
Gene_len =  upstream_3000bp.tolist() + exon1.tolist()+ intron1.tolist()+ exon2.tolist()+ intron2.tolist()+ exon3.tolist()+after_1000.tolist()+ before_1000.tolist()+exon_3.tolist()+ intron_2.tolist()+ exon_2.tolist()+ intron_1.tolist()+ exon_1.tolist()+downstream_3000.tolist()
print (Gene_len, len(Gene_len))
C_bw.close()

'''
x=range(1,170)
y=Gene_len.tolist()

fig = plt.figure()
plt.fill_between( x, y, color="skyblue", alpha=0.4)
plt.plot(x, y, color="Slateblue", alpha=0.6)
fig.savefig("foo.pdf", bbox_inches='tight')
'''
