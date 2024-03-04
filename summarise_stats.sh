while IFS="," read -r sample technology ignore
do
  echo -n "${sample},$(awk '$1=="Mean"&&($2=="coverage"||$2=="quality") {print $3;next} {if ($1=="N50") {print $2}}' alignments/${sample}.${technology}.mm2.stats | tr '\n' ',')"
  echo ""

done < <(tail +2 long_read_metadata.csv)
