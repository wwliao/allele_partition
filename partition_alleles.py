#!/usr/bin/env python3
import re
import argparse
from os.path import basename
from cyvcf2 import VCF, Writer

parser = argparse.ArgumentParser()
parser.add_argument("vcffile")
args = parser.parse_args()

invcf = VCF(args.vcffile)
prefix = re.search("(\S+)\.vcf\.gz", basename(args.vcffile))[1]
outvcf = Writer(f"{prefix}.partitioned.vcf.gz", invcf)
for record in invcf:
    clusters = []
    num_alleles = len(record.ALT)
    for count1, allele1 in sorted(zip(record.INFO.get("AC"), record.ALT), reverse=True):
        is_added = False
        # Why there are some alleles doesn't have any counts?
        if count1 > 0:
            for i, cluster in enumerate(clusters):
                for count2, allele2 in cluster:
                    if abs(len(allele2) - len(allele1)) < 50:
                        clusters[i].append((count1, allele1))
                        is_added = True
                        break
            if not is_added:
                clusters.append([(count1, allele1)])
   
    alleles = []
    counts = []
    an = record.INFO.get("AN")
    ns = record.INFO.get("NS")
    freqs = []
    for cluster in clusters:
        allele = cluster[0][1]
        count = sum(t[0] for t in cluster)
        alleles.append(allele)
        counts.append(count)
        freqs.append(f"{count/an:.6f}")
    alt = ",".join(alleles)
    ac = ",".join(map(str, counts))
    af = ",".join(freqs)
    record_str = f"{record.CHROM}\t{record.POS}\t{record.ID}\t{record.REF}\t{alt}\t60\t.\tAC={ac};AF={af};AN={an};NS={ns}"
    out_record = outvcf.variant_from_string(record_str)
    outvcf.write_record(out_record)
invcf.close()
outvcf.close()
