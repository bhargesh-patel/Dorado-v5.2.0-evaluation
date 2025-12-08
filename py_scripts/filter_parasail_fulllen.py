import pysam
import sys

input_file = sys.argv[1]
output_file = input_file.split(".sam")[0] + "_fulllen.sam"

samfile = pysam.AlignmentFile(input_file, "r")
outfile = pysam.AlignmentFile(output_file, 'w', template=samfile)

for read in samfile.fetch():
    if read.is_unmapped:
        continue

    ref_pos = read.reference_start
    cigar = read.cigartuples
    first_match_pos = None
    last_match_pos = ref_pos

    for op, length in cigar:
        if op in [0, 7, 8]:  # M, =, X
            if first_match_pos is None:
                first_match_pos = ref_pos
            ref_pos += length
            last_match_pos = ref_pos  # update last matched position
        elif op in [2, 3]:  # D, N
            ref_pos += length
        elif op in [1, 4, 5]:  # I, S, H
            continue  # doesn't advance ref
        else:
            continue

    reference_length = samfile.get_reference_length(read.reference_name)

    # Keep read if first match is near start and last match is near end
    if (first_match_pos is not None and first_match_pos < 25 and 
        (reference_length - last_match_pos) <= 20):
        outfile.write(read)

samfile.close()
outfile.close()
