#########################################################
# title:filter bam
# author:dawang
# date:20241023
#########################################################



import sys
import pysam
import pandas as pd
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
outsam = pysam.AlignmentFile(sys.argv[2], "wb",header=dict(samfile.header))

for i in samfile:
    try:
        if i.flag in [0,16]:
            if i.has_tag('GN'):
                outsam.write(i)
    except KeyError:
        continue

outsam.close()