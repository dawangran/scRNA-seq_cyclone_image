#########################################################
# title:chage bam tag
# author:dawang
# date:20230915
#########################################################



import sys
import pysam
samfile = pysam.AlignmentFile(sys.argv[1], "rb")
outsam = pysam.AlignmentFile(sys.argv[2], "wb",header=dict(samfile.header))
name = sys.argv[3]

for i in samfile:
    try:
        if i.has_tag('CB'):
            old = i.get_tag('CB')
            new = name + '_' + old
            i.set_tag('NB',new)
            outsam.write(i)
    except KeyError:
        continue

outsam.close()
