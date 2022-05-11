#!/usr/bin/env python3
import re
import argparse
from os.path import basename
from cyvcf2 import VCF, Writer

parser = argparse.ArgumentParser()
parser.add_argument("vcffile")
args = parser.parse_args()

invcf = VCF(args.vcffile)
invcf.add_info_to_header({"ID": "CLUSTER_AC", "Number": ".", "Type": "Integer", "Description": "Total number of alleles in called genotypes after clustering"})
invcf.add_info_to_header({"ID": "CLUSTER_AF", "Number": ".", "Type": "Float", "Description": "Estimated allele frequency in the range (0,1] after clustering"})
invcf.add_info_to_header({"ID": "CLUSTER_AI", "Number": ".", "Type": "String", "Description": "Representative allele index after clustering"})
prefix = re.search("(\S+)\.vcf\.gz", basename(args.vcffile))[1]
outvcf = Writer(f"{prefix}.partitioned.vcf.gz", invcf)
for record in invcf:
    clusters = []
    num_alt = len(record.ALT)
    if num_alt == 1:
        ref_count = record.INFO.get("AN") - record.INFO.get("AC")
        allele_list = sorted(zip([ref_count, record.INFO.get("AC")], range(num_alt + 1), [record.REF] + record.ALT), reverse=True)
    else:
        ref_count = record.INFO.get("AN") - sum(record.INFO.get("AC"))
        allele_list = sorted(zip([ref_count] + list(record.INFO.get("AC")), range(num_alt + 1), [record.REF] + record.ALT), reverse=True)
    for count1, idx1, allele1 in allele_list:
        is_added = False
        # Why there are some alleles doesn't have any counts?
        #if count1 > 0:
        for i in range(len(clusters)):
            for count2, idx2, allele2 in clusters[i]:
                if abs(len(allele2) - len(allele1)) < 50:
                    clusters[i].append((count1, idx1, allele1))
                    is_added = True
                    break
            if is_added:
                break
        if not is_added:
            clusters.append([(count1, idx1, allele1)])

    alleles = []
    counts = []
    freqs = []
    an = record.INFO.get("AN")
    for cluster in clusters:
        alleles.append("-".join([str(t[1]) for t in cluster]))
        count = sum(t[0] for t in cluster)
        counts.append(count)
        freqs.append(f"{count/an:.6f}")

    record.INFO["CLUSTER_AI"] = ",".join(alleles)
    record.INFO["CLUSTER_AC"] = ",".join(map(str, counts))
    record.INFO["CLUSTER_AF"] = ",".join(freqs)
    outvcf.write_record(record)
invcf.close()
outvcf.close()
