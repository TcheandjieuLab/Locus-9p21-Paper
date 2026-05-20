#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# User settings
# -----------------------------
PLINK2="plink2"

GENOTYPE_PREFIX="my_genotypes"
PHENO_FILE="phenotype.pheno"
PHENO_NAME="trait"
COVAR_FILE="covariates.cov"

CHR="9"
REGION_START="21000000"
REGION_END="23000000"

OUTPUT_PREFIX="conditional_analysis"
CONDITION_FILE="condition_list.txt"

N_RUNS=5

# Optional filters
MAF="0.05"
MAC="10"
MACH_R2="0.4"
MEMORY="5000"

# -----------------------------
# Initialize condition file
# -----------------------------
touch "$CONDITION_FILE"

# -----------------------------
# Iterative conditional analysis
# -----------------------------
for RUN in $(seq 1 "$N_RUNS"); do

    echo
    echo "========== RUN ${RUN} =========="

    RUN_OUTPUT="${OUTPUT_PREFIX}_run${RUN}"

    if [[ -s "$CONDITION_FILE" ]]; then
        CONDITION_ARG="--condition-list ${CONDITION_FILE}"
    else
        CONDITION_ARG=""
    fi

    $PLINK2 \
        --1 \
        --chr "$CHR" \
        --ci 0.95 \
        $CONDITION_ARG \
        --covar "$COVAR_FILE" \
        --covar-variance-standardize \
        --from-bp "$REGION_START" \
        --to-bp "$REGION_END" \
        --glm hide-covar cols=+a1freqcc,+ax,+a1freq,+machr2,+a1countcc firth-fallback \
        --mac "$MAC" \
        --mach-r2-filter "$MACH_R2" \
        --maf "$MAF" \
        --memory "$MEMORY" \
        --out "$RUN_OUTPUT" \
        --pfile "$GENOTYPE_PREFIX" \
        --pheno "$PHENO_FILE" \
        --pheno-name "$PHENO_NAME"

    RESULT_FILE=$(ls "${RUN_OUTPUT}".*.glm.* | head -n1)

    TOP_SNP=$(
        awk '
        NR==1{
            for(i=1;i<=NF;i++){
                if($i=="ID") id=i
                if($i=="P") p=i
                if($i=="TEST") test=i
            }
            next
        }
        ($p!="NA" && $p!="" && ($test=="" || $test=="ADD")) {
            print $id, $p
        }' "$RESULT_FILE" |
        sort -k2,2g |
        head -n1 |
        cut -d' ' -f1
    )

    if [[ -z "$TOP_SNP" ]]; then
        echo "No valid SNP found. Stopping."
        exit 0
    fi

    if grep -qx "$TOP_SNP" "$CONDITION_FILE"; then
        echo "Top SNP already conditioned on: $TOP_SNP"
        exit 0
    fi

    echo "$TOP_SNP" >> "$CONDITION_FILE"
    echo "Added top SNP: $TOP_SNP"

done

echo
echo "Finished ${N_RUNS} runs."
echo "Final conditioned variants:"
cat "$CONDITION_FILE"