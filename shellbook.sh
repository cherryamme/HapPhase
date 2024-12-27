### 
cargo build --release

/home/jiangchen/project/HapPhase/target/release/HapPhase --method cos -i "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam" -v "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/merge_output.vcf.gz" -r "chr5:70919543-70950000" -s 4

### 




### 
# Iter: i=05
workdir=/home/gushanshan/project/2024-11-27-SMN-FL-analysis
batch_num=115

sample=MCGD${batch_num}-$i
SMN1=$(/home/jiangchen/project/HapPhase/target/release/HapPhase -i $workdir/$sample/$sample.bam -v $workdir/$sample/merge_output.vcf.gz -r "chr5:70919543-70950000" -s 4|sed 's/Haplotype num: //g')
SMN2=$(/home/jiangchen/project/HapPhase/target/release/HapPhase -i $workdir/$sample/$sample.bam -v $workdir/$sample/merge_output.vcf.gz -r "chr5:70049543-70077704" -s 4|sed 's/Haplotype num: //g')
echo -e "$sample\t$SMN1\t$SMN2"
### 
