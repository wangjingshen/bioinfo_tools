import pysam
import sys

def reformat_bam(input_bam, output_bam):
    sam_in = pysam.AlignmentFile(input_bam, "rb")
    sam_out = pysam.AlignmentFile(output_bam, "wb", header=sam_in.header)

    for read in sam_in:
        if not (read.has_tag("CB") and read.has_tag("UB")):
            continue
        if not read.has_tag("NH"):
            continue
        nh = read.get_tag("NH")
        if nh != 1:
            continue
        cb = read.get_tag("CB")
        ub = read.get_tag("UB")
        old_qname = read.query_name
        new_qname = f"{cb}:{ub}:{old_qname}"
        read.query_name = new_qname
        sam_out.write(read)

    sam_in.close()
    sam_out.close()

if __name__ == "__main__":
    inbam = sys.argv[1]
    outbam = sys.argv[2]
    reformat_bam(inbam, outbam)