# -*- coding: utf-8 -*-
#!/usr/bin/env python3

import gzip
import sys
import time

start_time = time.time()

if len(sys.argv) != 3:
    sys.exit("Usage: python get_stat_fastq.py file1.fastq.gz file2.fastq.gz")

def process_file_pair(file1_path, file2_path):
    total_reads_1 = 0
    total_reads_2 = 0
    total_bases_1 = 0
    total_bases_2 = 0
    base_counts_1 = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
    base_counts_2 = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
    q20_bases_1 = 0
    q20_bases_2 = 0
    q30_bases_1 = 0
    q30_bases_2 = 0
    read_lengths_1 = []
    read_lengths_2 = []

    with gzip.open(file1_path, 'rb') as file1, gzip.open(file2_path, 'rb') as file2:
        while True:
            head1 = file1.readline().rstrip()
            head2 = file2.readline().rstrip()

            if not head1 or not head2:  
                break

            read1_head_num = head1.split()[1][0].startswith('1') # '1'
            read2_head_num = head2.split()[1][0].startswith('2') # '2'

            if read1_head_num and read2_head_num: 
                seq1 = file1.readline().rstrip()
                seq2 = file2.readline().rstrip()
                plus1 = file1.readline().rstrip()
                plus2 = file2.readline().rstrip()
                qual1 = file1.readline().rstrip()
                qual2 = file2.readline().rstrip()

                total_reads_1 += 1 
                read_lengths_1.append(len(seq1)) 
                total_reads_2 += 1 
                read_lengths_2.append(len(seq2)) 

                total_bases_1 += len(seq1) 
                total_bases_2 += len(seq2) 

                for q in qual1:
                    if ord(q) - 33 >= 20: # ord() : 
                        q20_bases_1 += 1
                    if ord(q) - 33 >= 30:
                        q30_bases_1 += 1

                for q in qual2:
                    if ord(q) - 33 >= 20:
                        q20_bases_2 += 1
                    if ord(q) - 33 >= 30:
                        q30_bases_2 += 1

                for base in seq1: 
                    if base in base_counts_1:
                        base_counts_1[base] += 1

                for base in seq2:
                    if base in base_counts_2:
                        base_counts_2[base] += 1
            else:
                if read1_head_num != '1' or read2_head_num != '2':
                    print("Read1과 Read2 는 서로 pair 관계가 아님")
                sys.exit(1)

    avg_read_length_1 = sum(read_lengths_1) / len(read_lengths_1)
    avg_read_length_2 = sum(read_lengths_2) / len(read_lengths_2)

    if total_reads_1 != total_reads_2:
        print("Read1과 Read2의 전체 read 개수가 일치하지 않음")
        sys.exit(1)

    return (
        total_reads_1, total_reads_2,
        total_bases_1, total_bases_2,
        base_counts_1, base_counts_2,
        q20_bases_1, q20_bases_2,
        q30_bases_1, q30_bases_2,
        avg_read_length_1, avg_read_length_2
        )

def print_summary(
        sample_name,
        total_reads_1, total_reads_2,
        total_bases_1, total_bases_2,
        base_counts_1, base_counts_2,
        q20_bases_1, q20_bases_2,
        q30_bases_1, q30_bases_2,
        avg_read_length_1, avg_read_length_2
        ):

    q_20p_1 = round(q20_bases_1 * 100.0 / total_bases_1, 2)
    q_20p_2 = round(q20_bases_2 * 100.0 / total_bases_2, 2)
    q_20p_total = round((q_20p_1 + q_20p_2)/2, 2)

    q_30p_1 = round(q30_bases_1 * 100.0 / total_bases_1, 2)
    q_30p_2 = round(q30_bases_2 * 100.0 / total_bases_2, 2)
    q_30p_total = round((q_30p_1 + q_30p_2)/2, 2)

    Np_1 = round(base_counts_1['N'] * 100.0 / total_bases_1, 2)
    Np_2 = round(base_counts_2['N'] * 100.0 / total_bases_2, 2)
    Np_total = round((Np_1 + Np_2)/2, 2)
    gc_content_1 = round((base_counts_1['G'] + base_counts_1['C']) * 100.0 / total_bases_1, 2)
    gc_content_2 = round((base_counts_2['G'] + base_counts_2['C']) * 100.0 / total_bases_2, 2)
    gc_content_total = round((gc_content_1 + gc_content_2)/2 ,2)

    print("{} / {:,} / {:,} / {} / {} / {} / {}".format(sample_name,
                                                        total_bases_1 + total_bases_2,
                                                        total_reads_1 + total_reads_2,
                                                        Np_total,
                                                        gc_content_total,
                                                        q_20p_total,
                                                        q_30p_total
                                                        ))

    print "SampleName :", sample_name
    for base in ['A', 'T', 'G', 'C', 'N']:
        print("{} base : {:,}".format(base, base_counts_1[base] + base_counts_2[base]))
    print "Q20 Bases :", "{:,}".format(q20_bases_1 + q20_bases_2)
    print "Q30 Bases :", "{:,}".format(q30_bases_1 + q30_bases_2)
    print("Avg Read Length: {:.2f}".format((avg_read_length_1 + avg_read_length_2) / 2))  
    print("------------------------------------------------------")

    print("{} / {:,} / {:,} / {} / {} / {} / {}".format(sample_name, total_bases_1, total_reads_1, Np_1, gc_content_1, q_20p_1, q_30p_1))
    print "SampleName :", sample_name + "_R1"
    for base in ['A', 'T', 'G', 'C', 'N']:
        print("{} base : {:,}".format(base, base_counts_1[base]))
    print "Q20 Bases :", "{:,}".format(q20_bases_1)
    print "Q30 Bases :", "{:,}".format(q30_bases_1)
    print("Avg Read Length: {:.2f}".format(avg_read_length_1))
    print("------------------------------------------------------")

    print("{} / {:,} / {:,} / {} / {} / {} / {}".format(sample_name, total_bases_2, total_reads_2, Np_2, gc_content_2, q_20p_2, q_30p_2))
    print "SampleName :", sample_name + "_R2"
    for base in ['A', 'T', 'G', 'C', 'N']:
        print("{} base : {:,}".format(base, base_counts_2[base]))
    print "Q20 Bases :", "{:,}".format(q20_bases_2)
    print "Q30 Bases :", "{:,}".format(q30_bases_2)
    print("Avg Read Length: {:.2f}".format(avg_read_length_2))
    print("------------------------------------------------------")

file_data = process_file_pair(sys.argv[1], sys.argv[2])

sample_name = sys.argv[1].split('_R1.fastq.gz')[0]

combined_total_reads_1 = file_data[0]  # Total Reads R1
combined_total_reads_2 = file_data[1]  # Total Reads R2
combined_total_bases_1 = file_data[2]  # Total Bases R1
combined_total_bases_2 = file_data[3]  # Total Bases R2
combined_base_counts_1 = file_data[4]  # Base Counts R1 ex) base_counts_1 = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
combined_base_counts_2 = file_data[5]  # Base Counts R2 ex) base_counts_2 = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
combined_q20_bases_1 = file_data[6]  # Q20 Bases R1
combined_q20_bases_2 = file_data[7]  # Q20 Bases R2
combined_q30_bases_1 = file_data[8]  # Q30 Bases R1
combined_q30_bases_2 = file_data[9]  # Q30 Bases R2
combined_avg_read_length_1 = file_data[10]  # Avg Read Length R1 ex) read_lengths_1 = []
combined_avg_read_length_2 = file_data[11]  # Avg Read Length R2 ex) read_lengths_2 = []

print_summary(
    sample_name,
    combined_total_reads_1, combined_total_reads_2,
    combined_total_bases_1, combined_total_bases_2,
    combined_base_counts_1, combined_base_counts_2,
    combined_q20_bases_1, combined_q20_bases_2,
    combined_q30_bases_1, combined_q30_bases_2,
    combined_avg_read_length_1, combined_avg_read_length_2
    )

end_time = time.time()
execution_time = end_time - start_time
print "걸린시간 :", execution_time, "초"
