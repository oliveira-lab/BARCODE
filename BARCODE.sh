#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob


#########################################################################
# File Name: BARCODE.sh
# Author(s): Angelina Beavogui and Vinicius S. Kavagutti
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 18/07/2025
#########################################################################


####### Initialization ####### 
eval "$(micromamba shell hook bash)"
micromamba activate BARCODE_pipeline

####### Fixed parameters #######
MIN_COMPL=90
MAX_CONT=5
MIN_GENOMES=${MIN_GENOMES:-2}
ANI_THRESH=1      
THREADS=${THREADS:-10}        

####### Usage & Options #######
usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT_DIR [-t THREADS] [-g MIN_GENOMES]

Description:
  This script processes genomic data through a pipeline involving quality control (CheckM),
  dereplication (dRep), annotation (Prokka), defense system analysis (DefenseFinder),
  and core/quasi-core analysis.

Options:
  -i  Path to the directory containing species folders (required)
  -t  Number of threads to use (default: $THREADS)
  -g  Minimum number of genomes to retain per species after QC and dereplication (default: $MIN_GENOMES)
  -h  Display this help message and exit

Example:
  $(basename "$0") -i /path/to/input -t 16 -g 3
EOF
  exit 1
}

INPUT_DIR=""
while getopts "i:t:g:h" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        g) MIN_GENOMES="$OPTARG" ;;
        *) usage ;;
    esac
done
[[ -d "$INPUT_DIR" ]] || { echo "Invalid INPUT_DIR"; usage; }


find "$INPUT_DIR" -type f -name "*.fasta" | while IFS= read -r genome_file; do
    new_file="${genome_file%.fasta}.fna"
    mv -- "$genome_file" "$new_file"
    echo "Renaming  ${genome_file##*/} → ${new_file##*/}"
done

for d in "$INPUT_DIR"/*/; do
  [[ -d $d ]] || continue
  base=${d%/}    
  base=${base##*/}

  if [[ "$base" == "${base//[^A-Za-z0-9]/__}" ]]; then
    continue
  fi

  new=$(printf '%s' "$base" | sed -E 's/[^A-Za-z0-9]+/__/g')
  new=$(printf '%s' "$new" | sed -E 's/(__)+/__/g')

  if [[ "$new" != "$base" ]]; then
    echo "Renaming $base → $new"
    mv -- "$d" "$INPUT_DIR/$new"
  fi
done


# Prepare output directory structure
RESULTS="$INPUT_DIR/results"
mkdir -p "$RESULTS"/{1_checkm,2_dRep,3_prokka,4_defense,rejected/{low_quality,too_few_genomes,dRep_too_few,duplicates}}


####### Quality Control with CheckM #######
run_checkm() {
  local species_dir=$1
  local sp=$(basename "$species_dir")
  local out="$RESULTS/1_checkm/$sp"
  mkdir -p "$out"

  echo "Running CheckM lineage_wf for $sp…"
  checkm taxonomy_wf domain Bacteria -x fna -t "$THREADS" "$species_dir" "$out"

  local stats="$out/storage/bin_stats_ext.tsv"
  [[ -s "$stats" ]] || { echo "$stats not found or empty"; return 1; }

  # Clear lists
  : > "$out/qualified.txt"
  : > "$out/rejected.txt"


  while IFS=$'\t' read -r id dict; do
    # Extract Completeness and Contamination
    c=$(grep -oP "'Completeness':\s*\K[0-9.]+" <<<"$dict" || echo 0)
    t=$(grep -oP "'Contamination':\s*\K[0-9.]+" <<<"$dict" || echo 999)

    # Check thresholds
    ok=$(awk -v c="$c" -v t="$t" -v MINC="$MIN_COMPL" -v MAXC="$MAX_CONT" \
               'BEGIN{print (c>=MINC && t<=MAXC)?1:0}')

    if [[ $ok -eq 1 ]]; then
      echo "$id" >> "$out/qualified.txt"
    else
      echo "$id" >> "$out/rejected.txt"
    fi
  done < "$stats"

  # Move rejected .fna files
  mkdir -p "$RESULTS/rejected/low_quality/$sp"
  while read -r g; do
    mv "$species_dir/$g.fna" "$RESULTS/rejected/low_quality/$sp"/
  done < "$out/rejected.txt" 2>/dev/null || true

  # Count qualified genomes
  mapfile -t passed < "$out/qualified.txt"
  local n=${#passed[@]}
  if (( n < MIN_GENOMES )); then
    echo "$sp has only $n QC-passed genomes (< $MIN_GENOMES) → rejecting species"
    mv "$species_dir" "$RESULTS/rejected/too_few_genomes/$sp"
    return 1
  fi

  return 0
}

####### Dereplication with dRep #######
run_drep() {
  local species_dir=$1
  local sp=$(basename "$species_dir")
  local out="$RESULTS/2_dRep/$sp"
  mkdir -p "$out"

  echo "Running dRep dereplicate for $sp…"
  dRep dereplicate "$out" -g "$species_dir"/*.fna -pa $ANI_THRESH --skip_plots --ignoreGenomeQuality

  local wdb="$out/data_tables/Wdb.csv"
  if [[ ! -f "$wdb" ]]; then
    echo "$wdb not found" >&2
    return 1
  fi

  mapfile -t keep < <(
    awk -F, 'NR>1 {
      fname = $1
      sub(/.*\//,"",fname)    
      sub(/\.fna$/,"",fname)  
      print fname
    }' "$wdb"
  )


  declare -A is_keep
  for id in "${keep[@]}"; do
    is_keep[$id]=1
  done

  # Move duplicates to rejected folder
  mkdir -p "$RESULTS/rejected/duplicates/$sp"
  for f in "$species_dir"/*.fna; do
    base=${f##*/}; base=${base%.fna}
    if [[ -z "${is_keep[$base]:-}" ]]; then
      mv "$f" "$RESULTS/rejected/duplicates/$sp"/
    fi
  done

  # Count remaining genomes
  shopt -s nullglob
  local files=( "$species_dir"/*.fna )
  local n=${#files[@]}
  if (( n < MIN_GENOMES )); then
    echo "After dRep, $sp has only $n genomes (< $MIN_GENOMES) → rejecting species"
    mv "$species_dir" "$RESULTS/rejected/dRep_too_few/$sp"
    return 1
  fi

  return 0
}


####### Annotation with Prokka #######
run_prokka() {
  local species_dir=$1
  local sp=$(basename "$species_dir")
  local out="$RESULTS/3_prokka/$sp"
  mkdir -p "$out"

  echo "Running Prokka for $sp…"
  for f in "$species_dir"/*.fna; do
    base=${f##*/}
    base=${base%.fna}
    prokka "$f" --outdir "$out/$base" --prefix "$base" --cpus "$THREADS"
  done
}

####### DefenseFinder #######
run_defense() {
  local prokka_dir=$1
  local sp=$(basename "$prokka_dir")
  local out="$RESULTS/4_defense/$sp"
  mkdir -p "$out"
  echo "Running DefenseFinder for $sp…"
  find "$prokka_dir" -name "*.faa" | while read -r faa; do
    sub=$(basename "$(dirname "$faa")")
    defense-finder run "$faa" -o "$out/$sub" --models-dir /opt/defense-finder-models
  done
}


####### Main processing loop #######
for species in "$INPUT_DIR"/*; do
    [[ -d "$species" ]] || continue
    [[ $(basename "$species") == "results" ]] && continue
    if run_checkm "$species"; then
        if run_drep "$species"; then
            run_prokka "$species"
            run_defense "$RESULTS/3_prokka/$(basename "$species")"
        fi
    fi
done

micromamba deactivate

####### Core/Quasi-core Analysis #######
micromamba activate BARCODE_analysis

 
parse_defense_systems(){
  local sys_file=$1
  perl - "$sys_file" <<'EOF'
use strict; use warnings;
my $f=shift;
open F,"<",$f or die $!;
<F>;  # skip header
while(<F>){
  chomp;
  my @c=split"\t";
  s/-/_/g for @c[0,1,2,4,5,6,7,8];
  my ($sys,$type,$sub,$beg,$end,$prot,$cnt,$name)=@c[0,1,2,4,5,6,7,8];
  my $pref=join"-",$sys,$type,$sub,$name;
  if($prot eq '1'){
    print "$pref $end\n";
  } else {
    for my $p (split/,/,$prot){
      $p=~s/^\s+|\s+$//g;
      print "$pref $p\n";
    }
  }
}
close F;
EOF
}

# Collect all Prokka .faa files into results/defense/*/

for species in "$RESULTS"/4_defense/*/; do
  sp=$(basename "$species")
  mkdir -p "$species"
  mapfile -t faas < <(find "$RESULTS"/3_prokka/"$sp"/ -type f -name '*.faa' 2>/dev/null)
  if (( ${#faas[@]} )); then
    echo " - $sp : copying ${#faas[@]} .faa files"
    cp "${faas[@]}" "$species"/
  else
    echo " No Prokka .faa files for $sp"
  fi
done


####### Parse defense systems #######

DEF_DIR="$INPUT_DIR/results/4_defense"
for species in "$DEF_DIR"/*/; do
  [[ -d $species ]] || continue
  sp=$(basename "$species")

  blast_dir="$species/blastall"
  pan="$blast_dir/pan_DSs-df.txt"
  mkdir -p "$blast_dir"; : > "$pan"
  
  # Replace hyphens with underscores in all .tsv files under this species
  for file in "$species"/*/*.tsv; do
    [[ -f $file ]] || continue
    sed 's/-/_/g' "$file" > tmpfile && mv tmpfile "$file"
  done
  

  for sysf in "$species"/*/*_defense_finder_systems.tsv; do
    [[ -f $sysf ]] || continue
    parse_defense_systems "$sysf" >> "$pan"
  done
done

# BLAST 
for species in "$RESULTS/4_defense"/*/; do
  [[ -d "$species" ]] || continue
  sp=$(basename "$species")
  blast_dir="$species/blastall"
  pan="$blast_dir/pan_DSs-df.txt"

  pushd "$blast_dir" >/dev/null
  sort -u pan_DSs-df.txt > 1 && mv 1 pan_DSs-df.txt

  awk '{print $2}' pan_DSs-df.txt | sort > tmp1
  cat ../*.faa > tmp2
  sed 's/-/_/g' tmp2 > tmpfile && mv tmpfile tmp2
  perl -ne 'if(/^>(\S+)/){$p=$i{$1}} $p?print:chomp; $i{$_}=1 if @ARGV' tmp1 tmp2 > tmp3
  makeblastdb -in tmp3 -parse_seqids -dbtype prot -out blastdb
  blastp -db blastdb -query tmp3 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out blast_results-df.txt -num_threads "$THREADS"
  awk 'BEGIN { FS = OFS = "\t"; cnt = 0 }{minlen = ($13 < $14 ? $13 : $14); cov    = ($4 / minlen) * 100; if ($3 >= 95 && cov >= 80 && $11 < 1e-4) { key = $1 OFS $2; keep[key] = 1; cnt++; q[cnt] = $1; s[cnt] = $2; record[cnt] = $0;}} END {for (i = 1; i <= cnt; i++) {fwd = q[i] OFS s[i]; rev = s[i] OFS q[i]; if (keep[fwd] && keep[rev]) { print record[i] }}}' blast_results-df.txt > tmp4
  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' tmp4 > blast_results-df.txt
  rm -fr blastdb.* tmp1 tmp2 tmp3

    : > g-and-g_list.txt
    for f in ../*.faa; do
      awk '{ if($0~">") print FILENAME, $1 }' "$f" >> g-and-g_list.txt
    done
    sed 's/-/_/g' g-and-g_list.txt > tmpfile && mv tmpfile g-and-g_list.txt
    sed 's/\.\.\///; s/\.faa//' g-and-g_list.txt |sed 's/>//g' | awk '{print $2,$1}' > 1 && mv 1 g-and-g_list.txt

    cat pan_DSs-df.txt | sort -k2 | cut -f2 -d' ' > tmp.list
    rg -j 10 -f tmp.list g-and-g_list.txt > g-and-g_list-df.txt
    sed -e 's/[ \t]\+/\t/g' pan_DSs-df.txt > tmpfile && mv tmpfile pan_DSs-df.txt
    sed -e 's/[ \t]\+/\t/g' g-and-g_list-df.txt > tmpfile && mv tmpfile g-and-g_list-df.txt
    sort -k1 g-and-g_list-df.txt > g-and-g_list-df.sorted.txt
    sort -k2 pan_DSs-df.txt | join -1 2 -2 1 - g-and-g_list-df.sorted.txt > pan_DSs-df.tmp.txt && mv pan_DSs-df.tmp.txt pan_DSs-df.txt
    rm -f g-and-g_list-df.sorted.txt
    rm -fr tmp.list


    python3 - "$PWD/blast_results-df.txt" "$PWD/clean_blast_results-df.txt" <<'EOF'
import sys, pandas as pd
infile, outfile = sys.argv[1], sys.argv[2]
df = pd.read_csv(infile, sep='\t', header=None,
  names=['query','subject','identity','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
df = df[df['query'] != df['subject']]
df['pair'] = df.apply(lambda r: ''.join(sorted([r['query'],r['subject']])), axis=1)
df = df.drop_duplicates('pair').drop(columns='pair')
df.to_csv(outfile, sep='\t', index=False, header=False)
EOF
  sort -k1 clean_blast_results-df.txt > 11 && mv 11 clean_blast_results-df.txt


  cat ../*/*_defense_finder_genes.tsv | cut -f2,3,13,18 | sed 's/ /\t/g' > cov.txt
  grep -v "hit_id" cov.txt | sort -k1 > tmp && mv tmp cov.txt
  

  cut -f1 cov.txt  | sort  | uniq -c | sort -n | awk '$1==2{print $2}' > dup.list

  if [[ -s dup.list ]]; then
    # Identify duplicates in cov.txt
    echo "$(wc -l < dup.list) duplicated IDs found"

    rg -j "$THREADS" -w -f dup.list clean_blast_results-df.txt > dup.blast || true

    awk '{print "s/" $1 "/" $1 "_dup/g"}' dup.list  | sed -f - dup.blast > dup.blast.renamed 

    cat clean_blast_results-df.txt dup.blast.renamed > clean_blast_results-df2.txt

    rm dup.blast dup.blast.renamed

  else
    cp clean_blast_results-df.txt clean_blast_results-df2.txt
  fi

  cat pan_DSs-df.txt | cut -f1,2 -d' ' > pan_DSs-df.tmp
  sed 's/[ \t]\+/\t/g' pan_DSs-df.tmp > tmpfile && mv tmpfile pan_DSs-df.tmp
  join -1 1 -2 1 -o 1.1,1.2,2.2 <(sort pan_DSs-df.tmp) <(sort -k1,1 g-and-g_list-df.txt) > merged_table.txt
  sed 's/[ \t]\+/\t/g' merged_table.txt > tmpfile && mv tmpfile merged_table.txt
  awk -F'\t' 'BEGIN{OFS="\t"} {if ($1 in seen) {seen[$1]++; $1 = $1 "_dup"} else seen[$1]++} 1' merged_table.txt > renamed_table.txt
  awk -F'\t' 'BEGIN{OFS="\t"} {if ($1 in seen) {seen[$1]++; $1 = $1 "_dup"} else seen[$1]++} 1' cov.txt > cov2.txt

    join -1 1 -2 1 -o 1.1,1.2,1.3,2.2,2.3,2.4 -t $'\t' renamed_table.txt cov2.txt > merged_table2.txt && \
    if [ $(wc -l < renamed_table.txt) -eq $(wc -l < cov2.txt) ]; then \
        echo "Quality check passed. Both files have the same number of lines."; \
    else \
        echo "Quality check failed. The number of lines in the files is different."; \
        echo "PWD: $(pwd)"; \
    fi


  rm -f dup.list dup.blast pan_DSs-df.tmp merged_table.txt renamed_table.txt cov2.txt 


  Rscript --vanilla - <<'EOF'
suppressPackageStartupMessages(library(dplyr))


df1 <- read.delim("clean_blast_results-df2.txt", header = FALSE, sep = "\t")
df2 <- read.delim("merged_table2.txt",          header = FALSE, sep = "\t")

res1 <- df1 %>% left_join(df2, by = c("V1" = "V1"))
res2 <- res1 %>% left_join(df2, by = c("V2.x" = "V1"))

# Write final table
write.table(res2,
            file      = "final_table-df.txt",
            sep       = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote     = FALSE)
EOF

  rm -f clean_blast* cov.txt g-and-g_list* merged_table2.txt
  popd >/dev/null
done

# Add species column
for species in "$RESULTS"/4_defense/*/; do
  sp=$(basename "$species")
  blast_dir="$species/blastall"
  infile="$blast_dir/final_table-df.txt"
  outfile="$blast_dir/${sp}_final_table-df.tsv"
  if [[ -f "$infile" ]]; then
    awk -v sp="$sp" -F'\t' 'BEGIN{OFS=FS}{
      for(i=1;i<=NF;i++){
        if($i~/^GCA_.*_genomic/){
          idx=index($i,"_genomic")
          pre=substr($i,1,idx-1); suf=substr($i,idx)
          gsub(/-/, "_", pre)
          $i=pre suf
        }
      }
      print $0, sp
    }' "$infile" > "$outfile"
  fi
done


####### Extract relevant columns and reformat #######
process_pidfiles() {
  for def_dir in "$RESULTS"/4_defense/*/; do
    sp=$(basename "$def_dir")
    blast_tsv="$def_dir/blastall/${sp}_final_table-df.tsv"
    core_dir="$RESULTS/5_core_qscore/$sp"

    mkdir -p "$core_dir"

      awk -F'\t' -v OFS='\t' '{print $1, $2, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23}' \
        "$blast_tsv" > "$core_dir/${sp}_final_table-df.3"

      awk -F'\t' '{
        print $1, $3, $4, $5, $6, $7, $13 > "temp1"
        print $2, $8, $9, $10, $11, $12, $13 > "temp2"
      }' "$core_dir/${sp}_final_table-df.3"
      cat temp1 temp2 > "$core_dir/${sp}_final_table-df.4"
      rm temp1 temp2
      
      awk -F'\t' -v OFS='\t' '{ gsub(/ +/,"\t"); print }' "$core_dir/${sp}_final_table-df.4" \
        | grep -v "accessory" > "$core_dir/${sp}_final_table-df.5"
        
      sort "$core_dir/${sp}_final_table-df.5" | uniq > "$core_dir/${sp}_final_table-df.6"
      
      awk -F'\t' -v OFS='\t' '{ gsub(/-/,"\t"); print }' "$core_dir/${sp}_final_table-df.6" > "$core_dir/${sp}_final_table-df.7"
      
      rm "$core_dir/${sp}_final_table-df.3" "$core_dir/${sp}_final_table-df.4" "$core_dir/${sp}_final_table-df.5" "$core_dir/${sp}_final_table-df.6"
    done
}


process_pidfiles

trim_duplicates() {
  local sp=$1
  local core_dir="$RESULTS/5_core_qscore/$sp"
  local def_dir="$RESULTS/4_defense/$sp"

  : > "$core_dir/tmp"
  : > "$core_dir/tmp1"

  for df in "$def_dir"/*/*_defense_finder_genes.tsv; do
    awk 'NR>1' "$df" >> "$core_dir/tmp"
  done

  sort -k2,2 -k17,17r "$core_dir/tmp" | sort -u -k2,2 | awk '{print $2,$3}' > "$core_dir/tmp1"
  sed 's/_dup//g' "$core_dir/${sp}_final_table-df.7" > tmpfile && mv tmpfile "$core_dir/${sp}_final_table-df.7"
  awk 'FNR==NR{a[$1]=$0;next}{print $0, a[$1]}' "$core_dir/tmp1" "$core_dir/${sp}_final_table-df.7" | awk '{if($7==$12) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' > "$core_dir/${sp}_final_table-df.7.trimmed"

}

for core_dir in "$RESULTS"/5_core_qscore/*/; do
  sp=$(basename "$core_dir")
  trim_duplicates "$sp"
done

####### Calculating presence/absence #######
run_fast_sys() {
  local core_dir=$1
  local sp=$(basename "$core_dir")
  local infile="$core_dir/${sp}_final_table-df.7.trimmed"
  local n_genomes
  # Count .fna in species folder
  n_genomes=$(ls -1 "$INPUT_DIR/$sp"/*.fna 2>/dev/null | wc -l)
  local cpus=$THREADS

  if [[ ! -f "$infile" ]]; then
    echo "No .7 files for $sp"
    return
  fi

  python3 - "$infile" "$n_genomes" "$cpus" <<'EOF'
import pandas as pd, sys
from multiprocessing import Pool

infile, n_genomes, cpus = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

df = pd.read_csv(infile, sep='\t', header=None, engine='python')


col6 = df[5].unique()
col4 = df[3].unique()

def proc(val4):
    return [ df[(df[3]==val4)&
                (df[5].str.contains(str(val6), na=False))].shape[0]
             for val6 in col6 ]

pool = Pool(cpus)
res = pool.map(proc, col4)
pool.close(); pool.join()

data = {'id': col6}
for v4, counts in zip(col4, res):
    data[v4] = counts
out = pd.DataFrame(data)

for c in out.columns[1:]:
    out[c] = (out[c]>0).astype(int)

sums = out.iloc[:,1:].sum()
tot = pd.DataFrame({'id':['TOTAL_COUNTS'], **{c:[sums[c]] for c in sums.index}})
out = pd.concat([out, tot], ignore_index=True)

if n_genomes>0:
    perc = {'id':['percentage']}
    for c in sums.index:
        perc[c] = [ round(sums[c]/n_genomes,3) ]
    out = pd.concat([out, pd.DataFrame(perc)], ignore_index=True)

out.to_csv(infile.rsplit('.',1)[0]+'.counts.tsv', sep='\t', index=False)
EOF
}

for core_dir in "$RESULTS"/5_core_qscore/*; do
  [[ -d "$core_dir" ]] || continue
  run_fast_sys "$core_dir"
done

for d in "$RESULTS"/4_defense/*/blastall; do
  rm -rf "$d"
done

####### Transpose counts tables #######
run_transpose() {
  local core_dir=$1
  local sp=$(basename "$core_dir")
  local infile="$core_dir/${sp}_final_table-df.7.counts.tsv"
  local outfile="$core_dir/${sp}_final_table-df.counts.transposed.tsv"

  if [[ ! -f "$infile" ]]; then
    echo " No *counts.tsv files for $sp"
    return
  fi

  python3 - "$infile" "$outfile" <<'EOF'
import sys

def transpose_columns(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    header = lines[0].rstrip('\n').split('\t')
    second_last = lines[-2].rstrip('\n').split('\t')
    last = lines[-1].rstrip('\n').split('\t')
    return list(zip(header, second_last, last))

def write_transposed(data, output_path):
    with open(output_path, 'w') as f:
        for row in data:
            f.write('\t'.join(row) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: transpose.py input_file output_file")
    inp, out = sys.argv[1], sys.argv[2]
    trans = transpose_columns(inp)
    write_transposed(trans, out)
EOF
}

for core_dir in "$RESULTS"/5_core_qscore/*/; do
  [[ -d "$core_dir" ]] || continue
  run_transpose "$core_dir"
done


####### Merge tables #######
run_merge() {
  local core_dir=$1
  local sp=$(basename "$core_dir")
  local file1="$core_dir/${sp}_final_table-df.7.trimmed"
  local file2="$core_dir/${sp}_final_table-df.counts.transposed.tsv"
  local outfile="$core_dir/${sp}_final_table-df.7.merge"

  if [[ ! -f "$file1" ]]; then
    echo " No .7 files for $sp"
    return
  fi
  if [[ ! -f "$file2" ]]; then
    echo " No counts.transposed.tsv files for $sp"
    return
  fi

  sed -e 's/\r$//' -e 's/\t$//' -e 's/\t\+/\t/g' "$file1" > tmpfile && mv tmpfile "$file1"

  python3 - "$file1" "$file2" "$outfile" <<'EOF'
import sys, pandas as pd

def merge_files(f1, f2, out):
    df2 = pd.read_csv(f2, sep='\t', header=0)
    mapping = df2.set_index('id')[['TOTAL_COUNTS','percentage']].to_dict(orient='index')

    df1 = pd.read_csv(f1, sep='\t', header=None)
    df1.columns = ["col1","col2","col3","col4","col5","col6","col7","col8","col9","col10"]

    df1['TOTAL_COUNTS'] = df1['col4'].map(lambda x: mapping.get(x,{}).get('TOTAL_COUNTS'))
    df1['percentage']   = df1['col4'].map(lambda x: mapping.get(x,{}).get('percentage'))

    df1.to_csv(out, sep='\t', index=False, header=False)

if __name__ == "__main__":
    if len(sys.argv)!=4:
        sys.exit("Usage: merge.py <in1> <in2> <out>")
    merge_files(*sys.argv[1:])
EOF
}

for core_dir in "$RESULTS"/5_core_qscore/*/; do
  [[ -d "$core_dir" ]] || continue
  run_merge "$core_dir"
done


####### Extract present defense system #######

for species_dir in "$RESULTS"/4_defense/*/; do
  sp=$(basename "$species_dir")
  core_dir="$RESULTS/5_core_qscore/$sp"
  n_genomes=$(ls -1 "$INPUT_DIR/$sp"/*.fna 2>/dev/null | wc -l)

  for f in "$species_dir"/*/*_defense_finder_systems.tsv; do
    [[ -f "$f" ]] || continue
    awk 'NR>1{print $7,$9}' "$f" >> "$core_dir/tmp2"
  done

  awk '{n=split($1,a,","); for(i=1;i<=n;i++) print a[i], "system" NR, $2}' "$core_dir/tmp2" > "$core_dir/tmp3"
  sort -k2,2 -k17,17r "$core_dir/tmp" | sort -u -k2,2 | awk '{print $2,$3,$13}' > "$core_dir/tmp4"
  awk 'FNR==NR{a[$1]=$0;next}{print $0, a[$1]}' "$core_dir/tmp4" "$core_dir/tmp3" | awk '$3 ~ $5' | awk '{if($0~"mandatory") print}' | awk '{count[$2]++; lines[NR]=$0; keys[NR]=$2} END {for (i=1; i<=NR; i++) print lines[i], count[keys[i]]}' | awk '{print $1,$2,$7}' > "$core_dir/tmp5"
  cut -f1,3,4,5,6,7,10,11,12 "$core_dir/${sp}_final_table-df.7.merge" > "$core_dir/${sp}_defense_systems_percentages_presence.tsv"
  awk 'FNR==NR{a[$1]=$0;next}{print $0, a[$1]}' "$core_dir/tmp5" "$core_dir/${sp}_defense_systems_percentages_presence.tsv" | awk '{if($9>=0.9) print}' | awk '{key = $11; genome[key] = $5; family[key] = $2; subfam[key] = $3; compo[key] = $4; species[key] = $7; number_df[key] = $12; ids[key] = (key in ids ? ids[key]","$1 : $1)} END {for (k in genome) {print genome[k]"\t"family[k]"\t"subfam[k]"\t"compo[k]"\t"species[k]"\t"ids[k]"\t"number_df[k]}}' | awk -F '\t' 'BEGIN{OFS="\t"} {split($6, a, ",");  count4 = length(a); print $0, count4}' > "$core_dir/tmp6"
  awk -F'\t' -v OFS='\t' \
      -v total_genomes="$n_genomes" -v sp="$sp" \
    '$NF == $(NF-1) {
       key = $3 "_" $1
       if (!(key in g_seen)) { g_seen[key]=1; count[$3]++ }
       lines[$3] = (lines[$3] ? lines[$3] RS $0 : $0)
     }
     END {
       for (s in lines) {
         if      (count[s] == total_genomes) filename = core_dir "/" sp "_core.txt"
         else if (count[s] >= 0.9*total_genomes) filename = core_dir "/" sp "_quasicore.txt"
         else continue
         print lines[s] > filename
       }
     }' core_dir="$core_dir" "$core_dir/tmp6"

  sed '1i Hit_ID\tType\tSubtype\tName_of_profiles_in_sys\tGenome_ID\tSpecies\tName_of_profiles_in_hit\tTotal_counts\tPercentage' "$core_dir/${sp}_defense_systems_percentages_presence.tsv" > tmpfile && mv tmpfile "$core_dir/${sp}_defense_systems_percentages_presence.tsv"
  rm -rf "$core_dir"/*.7* "$core_dir"/tmp* "$core_dir"/*transposed.tsv

  for suffix in core quasicore; do
    file="$core_dir/${sp}_${suffix}.txt"
    if [[ -f "$file" ]]; then
      cut -f1-6,9- "$file" > "$file.tmp" && mv "$file.tmp" "$file"
      sed '1i Genome_ID\tType\tSubtype\tName_of_profiles_in_sys\tSpecies\tHit_ID' "$file" > tmpfile && mv tmpfile "$file"
    fi
  done
done

echo; echo "Core/Quasi-core Analysis completed "

