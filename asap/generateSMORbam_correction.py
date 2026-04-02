#!/usr/bin/env python3
# encoding: utf-8
'''
asap.generateSMORbam -- Generate a corrected consensus BAM file for SMOR analysis.
Logic: Merges overlapping pairs, corrects mismatches based on Phred scores, 
and maintains alignment integrity via dynamic CIGAR strings.
'''

import sys
import os
import argparse
import logging
import pysam
from operator import attrgetter

def grouped_pairs(iterable):
    """Groups name-sorted reads into pairs."""
    it = iter(iterable)
    for x in it:
        try:
            yield x, next(it)
        except StopIteration:
            yield x, None

def _get_consensus(read, pair, fill_char):
    QUAL_DIFF_THRESHOLD = 10
    # Determine the total span of the fragment (Union)
    union_start = min(read.reference_start, pair.reference_start)
    union_end = max(read.reference_end, pair.reference_end)
    
    stats = {'total': 0, 'corrected': 0, 'ns': 0}

    def get_read_map(r):
        mapping = {}
        seq = r.query_sequence
        qual = r.query_qualities
        # Map reference position to (base, quality)
        for q_pos, r_pos in r.get_aligned_pairs():
            if r_pos is not None and q_pos is not None:
                mapping[r_pos] = (seq[q_pos], qual[q_pos])
        return mapping

    map1 = get_read_map(read)
    map2 = get_read_map(pair)

    consensus_seq = ""
    consensus_qual = []
    cigartuples = []
    
    # Track CIGAR operations to prevent coordinate shifting (False SNPs)
    curr_op = -1 # 0: Match/Mismatch (M), 2: Deletion (D)
    curr_len = 0

    def add_cigar(op):
        nonlocal curr_op, curr_len
        if op == curr_op:
            curr_len += 1
        else:
            if curr_op != -1:
                cigartuples.append((curr_op, curr_len))
            curr_op = op
            curr_len = 1

    # Iterate through the reference span
    for pos in range(union_start, union_end):
        b1, q1 = map1.get(pos, (None, None))
        b2, q2 = map2.get(pos, (None, None))

        # Scenario 1: Deletion in both reads relative to reference
        if b1 is None and b2 is None:
            add_cigar(2) # Record a Deletion (D)
            continue

        # Scenario 2: At least one read has a base (Match/Mismatch region)
        add_cigar(0) # Record a Match (M)
        
        if b1 and b2:
            stats['total'] += 1
            if b1 == b2:
                # Agreement: Sum qualities (cap at 60)
                consensus_seq += b1
                consensus_qual.append(min(q1 + q2, 60))
            else:
                # Mismatch: Apply Quality-based correction
                if q1 >= q2 + QUAL_DIFF_THRESHOLD:
                    consensus_seq += b1
                    consensus_qual.append(max(0, q1 - q2))
                    stats['corrected'] += 1
                elif q2 >= q1 + QUAL_DIFF_THRESHOLD:
                    consensus_seq += b2
                    consensus_qual.append(max(0, q2 - q1))
                    stats['corrected'] += 1
                else:
                    # Ambiguous: Mask with fill_char
                    consensus_seq += fill_char
                    consensus_qual.append(0)
                    stats['ns'] += 1
        elif b1:
            # Only Read 1 covers this tail
            consensus_seq += b1
            consensus_qual.append(q1)
        elif b2:
            # Only Read 2 covers this tail
            consensus_seq += b2
            consensus_qual.append(q2)

    # Finalize the CIGAR string
    if curr_len > 0:
        cigartuples.append((curr_op, curr_len))

    return (consensus_seq, consensus_qual, cigartuples, union_start, stats)

def _write_bam(samdata, out_file, fill_char):
    # Use a temporary name for sorting to avoid "file-in-use" indexing errors
    tmp_out = out_file + ".unsorted.tmp"
    outdata = pysam.AlignmentFile(tmp_out, "wb", template=samdata)
    
    grand_total_bases = 0
    grand_corrected = 0
    grand_ns = 0
    total_pairs = 0

    for ref_name in samdata.references:
        logging.info(f"Processing reference: {ref_name}")
        # Sorting by name is essential for find pairs in order
        reads = sorted(samdata.fetch(ref_name), key=attrgetter('query_name'))
        
        for read, pair in grouped_pairs(reads):
            if not pair or read.query_name != pair.query_name:
                continue
            if read.is_unmapped or pair.is_unmapped:
                continue
            
            # Simple overlap check
            if read.reference_end < pair.reference_start or pair.reference_end < read.reference_start:
                continue

            try:
                seq, qual, cigar, start, stats = _get_consensus(read, pair, fill_char)
                if seq:
                    new_read = pysam.AlignedSegment()
                    new_read.query_name = read.query_name
                    new_read.is_paired = False
                    new_read.reference_id = read.reference_id
                    new_read.reference_start = start
                    new_read.query_sequence = seq
                    new_read.query_qualities = qual
                    new_read.cigartuples = cigar
                    new_read.mapping_quality = max(read.mapping_quality, pair.mapping_quality)
                    outdata.write(new_read)
                    
                    grand_total_bases += stats['total']
                    grand_corrected += stats['corrected']
                    grand_ns += stats['ns']
                    total_pairs += 1
            except Exception as e:
                logging.error(f"Error processing {read.query_name}: {e}")

    outdata.close()
    
    # Summary Log
    logging.info("-" * 30)
    logging.info(f"Pairs Processed:  {total_pairs}")
    logging.info(f"Overlap Bases:    {grand_total_bases}")
    logging.info(f"Corrected:        {grand_corrected} ({(grand_corrected/max(1,grand_total_bases))*100:.2f}%)")
    logging.info(f"Masked (N):       {grand_ns} ({(grand_ns/max(1,grand_total_bases))*100:.2f}%)")
    logging.info("-" * 30)

    # Coordinate Sort and Index
    logging.info(f"Sorting and Indexing {out_file}...")
    pysam.sort("-o", out_file, tmp_out)
    pysam.index(out_file)
    
    if os.path.exists(tmp_out):
        os.remove(tmp_out)

def main():
    parser = argparse.ArgumentParser(description="SMOR Consensus Generator with Corrected CIGARs")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("-o", "--out", help="Output BAM file name")
    parser.add_argument("-c", "--fill-character", default="N", help="Character for ambiguous mismatches")
    # Added a logfile argument to match what Nextflow expects
    parser.add_argument("-l", "--logfile", default="smor_processing.log", help="Log file name")
    
    args = parser.parse_args()
    
    # Updated logging to write to BOTH the console (stderr) and a file
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s %(levelname)s: %(message)s',
        handlers=[
            logging.FileHandler(args.logfile),
            logging.StreamHandler(sys.stderr)
        ]
    )

    if not args.out:
        args.out = os.path.basename(args.bam).replace(".bam", "_SMOR.bam")

    with pysam.AlignmentFile(args.bam, "rb") as samdata:
        _write_bam(samdata, args.out, args.fill_character)
        
if __name__ == "__main__":
    main()