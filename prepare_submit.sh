
cd clair3
for F in *
do
  S=${F%_*}
  T=${F##*_}
  echo "moving ${S}_${T}/merge_output.gvcf.gz to ../to_transfer/${S}.${T}.g.vcf.gz"
  echo "md5sum ../to_transfer/${S}.${T}.g.vcf.gz > ../to_transfer/${S}.${T}.g.vcf.gz.md5"
done
cd ..

cd sniffles2
for F in *.snf
do
  echo "bgzip -@ 2 -c ${F} > ../to_transfer/${F}.gz"
  echo "md5sum ../to_transfer/${F}.gz > ../to_transfer/${F}.gz.md5"
done
cd ..
