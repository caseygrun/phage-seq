
import pandas as pd


# for each sample

# decide based on metadata whether to cut with alpaca or synthetic primers

# mkdir -p "$OUT_DIR_CUT/"{forward,reverse,info}
# mkdir -p "$OUT_DIR_CUT1/"{forward,reverse,info}
# mkdir -p "$OUT_DIR_CUT2/"{forward,reverse,info}

# Process Kaz lab samples
for read in "1n" "1p" "2n" "2p" "3n" "3p" "4n" "4p"; do
    PATH_FW="$DATA_PATH/forward/$read.fastq.gz"
    PATH_RV="$DATA_PATH/reverse/$read.fastq.gz"
    echo $PATH_FW
    echo $PATH_RV

    echo "Cutadapt 1..."

    cutadapt -g 'CGCTCAGGTTGCAGCTCGTGGAGTC' -G 'ATACGGCACCGGCGCACCACTAG' \
        --cores=0 \
        --discard-untrimmed \
        -o "$OUT_DIR_CUT1"/forward/"$read.fastq.gz" -p "$OUT_DIR_CUT1"/reverse/"$read.fastq.gz" \
        $PATH_FW $PATH_RV > "$OUT_DIR_CUT1/info/$read.txt"

    echo "Cutadapt 2..."
    cutadapt -g 'CGCTCAGGTTGCAGCTCGTGGAGTC' -G 'ATACGGCACCGGCGCACCACTAG' \
        --cores=0 \
        --discard-untrimmed \
        -o "$OUT_DIR_CUT2"/forward/"$read.fastq.gz" -p "$OUT_DIR_CUT2"/reverse/"$read.fastq.gz" \
        $PATH_RV $PATH_FW > "$OUT_DIR_CUT2/info/$read.txt"

    echo "cat forward..."
    cat "$OUT_DIR_CUT1/forward/$read.fastq.gz" "$OUT_DIR_CUT2/forward/$read.fastq.gz" > \
        "$OUT_DIR_CUT/forward/$read.fastq.gz"

    echo "cat reverse..."
    cat "$OUT_DIR_CUT1/reverse/$read.fastq.gz" "$OUT_DIR_CUT2/reverse/$read.fastq.gz" > \
        "$OUT_DIR_CUT/reverse/$read.fastq.gz"
done

# Process Ring lab samples
for read in "5n" "5p" "6n" "6p" "7n" "7p" "8n" "8p"; do
    PATH_FW="$DATA_PATH/forward/$read.fastq.gz"
    PATH_RV="$DATA_PATH/reverse/$read.fastq.gz"
    echo $PATH_FW
    echo $PATH_RV

    echo "Cutadapt 1..."

    cutadapt -g 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNnagctgcgcggcgagc' \
        -G 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNCGTATGGGTAagaaccaccgct' \
        --cores=0 \
        --discard-untrimmed \
        -o "$OUT_DIR_CUT1"/forward/"$read.fastq.gz" -p "$OUT_DIR_CUT1"/reverse/"$read.fastq.gz" \
        $PATH_FW $PATH_RV > "$OUT_DIR_CUT1/info/$read.txt"

    echo "Cutadapt 2..."
    cutadapt -g 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNnagctgcgcggcgagc' \
        -G 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGNNNNNNNNCGTATGGGTAagaaccaccgct' \
        --cores=0 \
        --discard-untrimmed \
        -o "$OUT_DIR_CUT2"/forward/"$read.fastq.gz" -p "$OUT_DIR_CUT2"/reverse/"$read.fastq.gz" \
        $PATH_RV $PATH_FW > "$OUT_DIR_CUT2/info/$read.txt"

    echo "cat forward..."
    cat "$OUT_DIR_CUT1/forward/$read.fastq.gz" "$OUT_DIR_CUT2/forward/$read.fastq.gz" > \
        "$OUT_DIR_CUT/forward/$read.fastq.gz"

    echo "cat reverse..."
    cat "$OUT_DIR_CUT1/reverse/$read.fastq.gz" "$OUT_DIR_CUT2/reverse/$read.fastq.gz" > \
        "$OUT_DIR_CUT/reverse/$read.fastq.gz"
done

} >>"$LOG_PATH" 2>&1
