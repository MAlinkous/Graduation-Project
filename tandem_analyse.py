#!/usr/bin/python

import os
import re
from xml.dom.minidom import parse
import xml.dom.minidom
import operator
root = os.getcwd()
list=[]
for i in os.listdir(root):
    if os.path.splitext(i)[1] == ".xml":
        list.append(i)
for i in list:
    tree = xml.dom.minidom.parse(i)
    root = tree.documentElement
    source=root.getAttribute('label')
    label = re.sub(r'models from ', "", source)
    label=re.sub('\'','',label)
    p = re.compile(r':reversed$')
    fo=open('analyse.xlsx','w')
    fo.write('group_id')
    fo.write('\t')
    fo.write('group_expect')
    fo.write('\t')
    fo.write('protein_expect')
    fo.write('\t')
    fo.write('protein_label')
    fo.write('\t')
    fo.write('protein_sumI')
    fo.write('\t')
    fo.write('peptide_sumI')
    fo.write('\t')
    fo.write('hyperscore')
    fo.write('\t')
    fo.write('seq')
    fo.write('\t')
    fo.write("is_reversed")
    fo.write('\n')
    k=0
    for node in root.childNodes:
        if node.nodeType != 3 and node.getAttribute('type')=='model':
            pep_sumI=node.getAttribute('sumI')
            group_id=node.getAttribute('id')
            group_expect=node.getAttribute('expect')
            protein = node.getElementsByTagName("protein")
            for pro in protein:
                fo.write(group_id)
                fo.write('\t')
                fo.write(group_expect)
                fo.write('\t')
                fo.write(pro.getAttribute('expect'))
                fo.write('\t')
                fo.write(pro.getAttribute('label'))
                fo.write('\t')
                fo.write(pro.getAttribute('sumI'))
                fo.write('\t')
                fo.write(pep_sumI)
                fo.write('\t')
                peptide=pro.getElementsByTagName('domain')
                fo.write(peptide[0].getAttribute('hyperscore'))
                fo.write('\t')
                fo.write(peptide[0].getAttribute('seq'))
                fo.write('\t')
                note=pro.getElementsByTagName('note')
                text=note[0].firstChild.data
                if p.search(text)==None:
                    fo.write('0')
                else:
                    fo.write('1')
                fo.write('\n')
    fo.close()

    fi2=open('analyse.xlsx','r')
    fo2=open('res.txt','w')
    lines1=fi2.readlines()
    dic={}
    dic2={}
    for i in range(1,len(lines1)):
        s1=lines1[i].strip().split('\t')
        m=s1[0]
        dic2[i]=m
        dic[m]=s1[6]
    sorted_dic=sorted(dic.items(),key=operator.itemgetter(1),reverse=True) 
    fo2.write(lines1[0])
    for i in sorted_dic:
        for j in dic2:
            if i[0]==dic2[j]:
                fo2.write(lines1[j])
                break
    fi2.close()
    fo2.close()
    filename_in=label+"_output.txt"
    fi3=open('res.txt','r')
    fo3=open(filename_in,'w')
    lines2=fi3.readlines()
    fo3.write(lines1[0])
    def get_i():
        j=0.0
        for i in range(1,len(lines2)):
            s2=lines2[i].strip().split('\t')
            if s2[8]=='1':
                j=j+1
            fdr=(2*j)/i
            if fdr>=0.01:
                x=float(s2[6])
                return x
                break
    score=get_i()
    for i in lines1:
        s3=i.strip().split('\t')
        if s3[0]!='group_id' and float(s3[6])>score:
            fo3.write(i)
    fi3.close()
    fo3.close()
    print(filename_in)
