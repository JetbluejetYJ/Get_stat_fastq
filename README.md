
# 📊 Get_stat_fastq

## 📝 Overview

`get_stat_fastq.py` is a Python script for rapid quality control (QC) statistics of paired-end NGS data in compressed FASTQ format (`.fastq.gz`).  
It processes two paired FASTQ files (R1 and R2) and outputs essential sequencing statistics:

- **Total number of reads and bases**
- **Base composition** (A, T, G, C, N)
- **Q20 and Q30 base percentages**
- **N base percentage**
- **GC content percentage**
- **Average read length**

---

## 🚀 Usage

```bash
python get_stat_fastq.py sample_R1.fastq.gz sample_R2.fastq.gz
```
- **Input:** Two compressed FASTQ files (R1 and R2)
- **Output:** QC statistics printed to the terminal
---

## 🖨️ Output Example

```
Sample1 / 1,234,567,890 / 12,345,678 / 0.01 / 48.23 / 98.76 / 95.43
SampleName : Sample1
A base : 300,000,000
T base : 300,000,000
G base : 317,283,945
C base : 317,283,945
N base : 0
Q20 Bases : 1,219,876,543
Q30 Bases : 1,176,543,210
Avg Read Length: 100.00
------------------------------------------------------
Sample1 / 617,283,945 / 6,172,839 / 0.01 / 48.23 / 98.76 / 95.43
SampleName : Sample1_R1
A base : 150,000,000
T base : 150,000,000
G base : 158,641,972
C base : 158,641,973
N base : 0
Q20 Bases : 609,938,271
Q30 Bases : 588,271,605
Avg Read Length: 100.00
------------------------------------------------------
Sample1 / 617,283,945 / 6,172,839 / 0.01 / 48.23 / 98.76 / 95.43
SampleName : Sample1_R2
A base : 150,000,000
T base : 150,000,000
G base : 158,641,973
C base : 158,641,972
N base : 0
Q20 Bases : 609,938,272
Q30 Bases : 588,271,605
Avg Read Length: 100.00
------------------------------------------------------
Elapsed time : 12.34 seconds
```

---

## 🧮 Calculation Methods

### 1. Total Reads and Bases

- **Total Reads:** Number of read records in each FASTQ file.
- **Total Bases:** Sum of the lengths of all sequence lines.

### 2. Base Composition

For each base (A, T, G, C, N):

$$
\text{Base Count} = \sum_{i=1}^{N} \text{count of base in read}_i
$$

### 3. Q20 and Q30 Bases

- **Phred Score Calculation:**  
  $Q = \text{ord}(q) - 33$  
  where $q$ is the ASCII character in the quality string.

- **Q20/Q30 Percentage:**

$$
\text{Q20\%} = \frac{\text{Number of Q} \geq 20 \text{ bases}}{\text{Total bases}} \times 100
$$

$$\text{Q30\%} = \frac{\text{Number of Q} \geq 30 \text{ bases}}{\text{Total bases}} \times 100$$

### 4. N Base Percentage

$$
\text{N\%} = \frac{\text{Number of N bases}}{\text{Total bases}} \times 100
$$

### 5. GC Content

$$
\text{GC\%} = \frac{\text{Number of G bases} + \text{Number of C bases}}{\text{Total bases}} \times 100
$$

### 6. Average Read Length

$$
\text{Average Read Length} = \frac{\sum_{i=1}^{N} \text{Length of read}_i}{N}
$$

---

## 🛠️ Script Workflow

1. **Input Validation:** Checks that two input files are provided.
2. **File Processing:** Reads both files in parallel, line by line, ensuring paired reads.
3. **Statistics Calculation:** For each read, updates counts for bases, Q20/Q30, and read lengths.
4. **Summary Output:** Prints overall and per-file statistics in a human-readable format.
5. **Error Handling:** Exits if files are not properly paired or if read counts do not match.

---

## ⚠️ Notes

- Assumes standard Illumina FASTQ format and Phred+33 encoding.
- Both input files must be properly paired and gzip-compressed.
- Intended for quick QC, not for in-depth quality analysis.

---
