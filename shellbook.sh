### 
cargo build --release

/home/jiangchen/project/HapPhase/target/release/HapPhase --method cos -i "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam" -v "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/merge_output.vcf.gz" -r "chr5:70919543-70950000" -s 4

### 

### 


bcftools view -T ^repeats.bed -o filtered_output.vcf -O v /home/gushanshan/project/2024-11-27-SMN-FL-analysis/MCGD115-06/merge_output.vcf.gz
bedtools intersect -a SMN.bed -b repeats.bed > SMN_filter.bed


### 


### 
# Iter: i=05
# Iter: i=05
for i in 1;do
workdir=/home/gushanshan/project/2024-11-27-SMN-FL-analysis
batch_num=102

sample=MCGD${batch_num}-$i
SMN1=$(/home/jiangchen/project/HapPhase/target/release/HapPhase -i $workdir/$sample/$sample.bam -v $workdir/$sample/merge_output.vcf.gz -b "/home/jiangchen/project/HapPhase/SMN_filter.bed" -r "chr5:70919543-70950000" -s 5|sed 's/Haplotype num: //g'|sed 's/Record num://g')
SMN2=$(/home/jiangchen/project/HapPhase/target/release/HapPhase -i $workdir/$sample/$sample.bam -v $workdir/$sample/merge_output.vcf.gz -b "/home/jiangchen/project/HapPhase/SMN_filter.bed" -r "chr5:70049543-70077704" -s 5|sed 's/Haplotype num: //g'|sed 's/Record num://g')

# 使用 IFS 和数组进行分割
IFS=$'\t' read -r -a SMN1_array <<< "$SMN1"
IFS=$'\t' read -r -a SMN2_array <<< "$SMN2"

# 获取record数量大小
SMN1_hapnum=${SMN1_array[0]:0:1}
SMN2_hapnum=${SMN2_array[0]:0:1}
SMN1_count=${SMN1_array[1]}
SMN2_count=${SMN2_array[1]}

# echo -e "raw: $sample\t$SMN1\t$SMN2"
# 比较并输出结果
if (( SMN1_count >= 5 * SMN2_count )); then
    echo -e "$sample\t${SMN1_array[0]}\t0"
elif (( SMN1_count <= SMN2_count / 5 )); then
    echo -e "$sample\t0\t$SMN2_hapnum"
elif (( SMN1_hapnum > 2 * SMN2_hapnum && SMN1_count <= 2 * SMN2_count )); then
	echo -e "$sample\t${SMN1_array[0]} (may error)\t${SMN2_array[0]}"
elif (( SMN2_hapnum > 2 * SMN1_hapnum && SMN2_count <= 2 * SMN1_count )); then
	echo -e "$sample\t${SMN1_array[0]}\t${SMN2_array[0]} (may error)"
else
    echo -e "$sample\t${SMN1_array[0]}\t${SMN2_array[0]}"
fi

done | code -
### 


### 

/home/jiangchen/project/HapPhase/target/release/HapPhase -i /home/gushanshan/project/2024-11-27-SMN-FL-analysis/MCGD115-06/MCGD115-06.bam -v /home/jiangchen/project/HapPhase/filtered_output.vcf.gz -r "chr5:70919543-70950000" -s 4|sed 's/Haplotype num: //g'
/home/jiangchen/project/HapPhase/target/release/HapPhase -i /home/gushanshan/project/2024-11-27-SMN-FL-analysis/MCGD115-06/MCGD115-06.bam -v /home/jiangchen/project/HapPhase/filtered_output.vcf.gz -r "chr5:70049543-70077704" -s 4|sed 's/Haplotype num: //g'
/home/jiangchen/project/HapPhase/target/release/HapPhase -i /home/gushanshan/project/2024-11-27-SMN-FL-analysis/MCGD115-06/MCGD115-06.bam -v /home/gushanshan/project/2024-11-27-SMN-FL-analysis/MCGD115-06/merge_output.vcf.gz -r "chr5:70919543-70950000" -s 4|sed 's/Haplotype num: //g'

### 
