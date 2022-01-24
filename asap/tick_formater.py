
#!/usr/bin/env python3
# encoding: utf-8
from lxml import etree as ET
import copy
import sys
import os
import argparse
from collections import defaultdict
from collections import Counter


def parse(read_dir, assay_name, sample_names):
    assay_file = "../fasta/" + assay_name + ".fasta"
    a = open(assay_file, 'a+')
    toSort = []
    for sample_name in sample_names:
        reads_file = read_dir + '/' + sample_name + '_' + assay_name + '.fasta'
        f = open(reads_file, 'r')
        c = Counter()
        total = 0
        for line in f:
            line = line.strip()
            if line[0] == '>':
                continue
            else:
                c[line] += 1
                total += 1
        passed = 0
        for allele in c.most_common():
            if passed >= 2 and (allele[1]/total) < .02:
                break
            p = []
            sample_name_trim = sample_name.split('-bowtie')[0]
            p.append(">"+sample_name_trim+"_"+str(allele[1])+"\n")
            p.append(allele[0] + "\n")
            p.append(allele[1])
            toSort.append(p)
            passed += 1
    toSort.sort(key=lambda x: int(x[2]), reverse=True)
    for allele in toSort:
        a.write(allele[0])
        a.write(allele[1])
    return

def main(argv=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-x", "--xml", help="Already reformatted ASAP output XML file. ")
    parser.add_argument("-o","--output", help="file to write new XML to. ")
    parser.add_argument("--switch", action='store_true', help="make parse fasta instead of create xml")
    parser.add_argument("-q", "--dir")
    parser.add_argument("-a", "--assay")
    parser.add_argument("-s", "--sample", nargs='*')
    args = parser.parse_args()
    if args.switch:
        read_dir = args.dir
        assay_name = args.assay
        sample_names = args.sample
        parse(read_dir, assay_name, sample_names)
        return
    else:
        out_file = args.output
        tree = ET.parse(args.xml)
        root = tree.getroot()

    #get the sequences, counts, and sample counts for each assay for assay-specific output
    all_dicts = defaultdict(defaultdict)
    for assay in root:
        if assay.get('type') == 'mixed' or assay.get('type') == 'ROI':
            currentAssayDict = defaultdict(list)
            for sample in assay:
                num_alleles = 0
                for allele in sample.iter('allele_sequence'):
                    if currentAssayDict[allele.text] == []:
                        li = [0, 0]
                        li[0] += int(allele.get('count'))
                        li[1] += 1
                        li.append([sample.get('name'), allele.get('count'), allele.get('percent')])
                        currentAssayDict[allele.text] = li
                    else:
                        lis = currentAssayDict[allele.text]
                        lis[0] += int(allele.get('count'))
                        lis[1] += 1
                        lis.append([sample.get('name'), allele.get('count'), allele.get('percent')])
                        currentAssayDict[allele.text] = lis
                    num_alleles += 1
                if num_alleles == 0: #to check for case where indel caused no 'allele_sequence's, even though reads aligned
                    if sample.find('./amplicon').get('reads') != "0":
                        consensus_string = "Something went wrong with getting consensus allele"
                        for consensus in sample.findall('./amplicon/consensus_sequence'):#there will only be one, but still easiest to loop
                            consensus_string = consensus.text
                        if currentAssayDict[consensus_string] == []:
                            li = [0, 0]
                            li[0] += int(sample.find('./amplicon').get('reads'))
                            li[1] += 1
                            li.append([sample.get('name') + '_CONSENSUS', sample.find('./amplicon').get('reads'), "100.0"])
                            currentAssayDict[consensus_string] = li
                        else:
                            lis = currentAssayDict[consensus_string]
                            lis[0] += int(sample.find('./amplicon').get('reads'))
                            lis[1] += 1
                            lis.append([sample.get('name') + '_CONSENSUS', sample.find('./amplicon').get('reads'), "100.0"])
                            currentAssayDict[consensus_string] = lis
                        if consensus_string == "Something went wrong with getting consensus allele": #double checking for case where there is no consensus to be found, this sample should be skipped because did not have reads aligning to assay
                            continue
            all_dicts[assay.get('name')] = currentAssayDict

    #write out this info to xml tree
    assay_counts_node = ET.SubElement(root, 'assay_counts')
    for dict in all_dicts.items():
        assay_node = ET.SubElement(assay_counts_node, 'assay')
        assay_node.set('name', dict[0])
        for allele in dict[1]:
            allele_node = ET.SubElement(assay_node, 'allele')
            allele_node.set('sequence',allele)
            allele_node.set('read_count', str(dict[1][allele][0]))
            allele_node.set('sample_count', str(dict[1][allele][1]))
            allele_node.set('hash', str(hash(allele)))
            for li in range(2, len(dict[1][allele])):
                sample_node = ET.SubElement(allele_node, "sample")
                sample_node.set('name', dict[1][allele][li][0])
                sample_node.set('count', dict[1][allele][li][1])
                sample_node.set('percent', dict[1][allele][li][2])


    #create a portion of xml for each allele's info
    for dict in all_dicts.items():
        for allele in dict[1]:
            allele_node = ET.SubElement(root, "allele_occurences")
            allele_node.set('sequence', allele)
            allele_node.set('hash', str(hash(allele)))
            allele_node.set('assay', str(dict[0]))
            for li in range(2, len(dict[1][allele])):
                sample_node = ET.SubElement(allele_node, "sample")
                sample_node.set('name', dict[1][allele][li][0])
                sample_node.set('count', dict[1][allele][li][1])
                sample_node.set('percent', dict[1][allele][li][2])

    tree.write(out_file, pretty_print=True)


if __name__ == "__main__":
    main()
