use std::collections::HashMap;
use log::info;
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::ext::BamRecordExtensions;





#[derive(Debug, Clone)]
pub struct HapArray {
    pub array: Vec<u8>,
    pub array_proportion: f64,
    pub count: usize,
}



#[derive(Debug, Clone)]
pub struct HapArrays {
    snp: Vec<(String, usize, char, char)>,
    array: Vec<HapArray>,
    length: usize,
    total_count: usize,

}
// 迭代这个结构时迭代内部的array字段
impl IntoIterator for HapArrays {
    type Item = HapArray;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.array.into_iter()
    }
}







pub fn get_bam_dict(input_bam: &str, snps: &Vec<(String, i64, char, char, f64, f64)>) -> HashMap<String, HashMap<(String, usize, char, char), bool>> {
    let start_time = std::time::Instant::now();
    let mut bam: bam::Reader = bam::Reader::from_path(input_bam).expect("Error opening BAM file");
    let header = bam.header().to_owned();
    let mut read_snp_map: HashMap<String, HashMap<(String, usize, char, char), bool>> = HashMap::new();
    let mut record_count = 0;
    let mut map_count = 0;
    // info!("SNPs: {:?}", snps);
    for record in bam.records() {
        let record = record.expect("Error reading BAM record");
        record_count += 1;
        if record.is_unmapped() {
            continue;
        }
        map_count += 1;
        let read_name = String::from_utf8(record.qname().to_vec()).expect("Error converting read name");
        let chrom = header.tid2name(record.tid() as u32);
        let chrom_str = String::from_utf8(chrom.to_vec()).expect("Error converting chromosome name");
        // TODO 使用rust_htslib 的read_pos方法进行优化
        // TODO 支持输入indel类型的snp
        for snp in snps {
            let (snp_chrom, snp_pos, snp_ref, snp_alt, _, _) = snp;
            if chrom_str == *snp_chrom {
                if let Ok(Some(read_pos)) = record.cigar().read_pos(*snp_pos as u32, true, false) {
                    let read_pos_usize = read_pos as usize - 1;
                    let read_base = record.seq().as_bytes()[read_pos_usize] as char;
                    let is_snp_present = read_base == *snp_alt;
                    read_snp_map.entry(read_name.clone())
                    .or_insert_with(HashMap::new)
                    .insert((snp_chrom.to_string(), *snp_pos as usize, *snp_ref, *snp_alt), is_snp_present);
                    // NOTE 用于 debug
                    // if *snp_pos == 70935946 {
                    //     info!("SNP: {:?}", snp);
                    //     info!("chrom:{:?} read_base:{:?}", snp_chrom.to_string(),  String::from_utf8_lossy(&record.seq().as_bytes()[read_pos_usize - 1..read_pos_usize + 2]));
                    //     info!("chrom:{:?} read_base:{:?}", snp_chrom.to_string(),  record.seq().as_bytes()[read_pos_usize - 1] as char);
                    // }
                }
            }
        }

    }
    // 删除没有覆盖到那些snp的reads
    let snp_count = snps.len();
    read_snp_map.retain(|_, snp_map| snp_map.len() >= snp_count);
    let elapsed_time = start_time.elapsed();
    info!("Get_bam_dict Elapsed time: {:?}", elapsed_time);
    info!("Total records: {}, map_count: {}, contain SNPs count: {}", record_count, map_count, read_snp_map.len());
    read_snp_map
}



use std::collections::{HashSet, VecDeque};
use std::collections::hash_map::Entry;
use std::collections::BTreeMap;
use std::cmp::Ordering;

pub fn process_mutation_data(read_mutation_presence: &HashMap<String, HashMap<(String, usize, char, char), bool>>) -> HapArrays {
    // let start_time = std::time::Instant::now();
    // Collect all SNP keys and sort them to ensure consistent ordering
    let mut snp_keys: Vec<(String, usize, char, char)> = read_mutation_presence
        .values()
        .flat_map(|snp_map| snp_map.keys().cloned())
        .collect();
    snp_keys.sort();
    snp_keys.dedup();

    // Convert bool to u8 and form vectors
    let mut converted_data: HashMap<String, Vec<u8>> = HashMap::new();
    for (read_name, snp_map) in read_mutation_presence {
        let values: Vec<u8> = snp_keys
            .iter()
            .map(|key| if *snp_map.get(key).unwrap_or(&false) { 1 } else { 0 })
            .collect();
        converted_data.insert(read_name.clone(), values);
    }

    // Count occurrences of each vector
    let mut array_count: HashMap<Vec<u8>, usize> = HashMap::new();
    for array in converted_data.values() {
        *array_count.entry(array.clone()).or_insert(0) += 1;
    }

    // Calculate the total count
    let total_count: usize = array_count.values().sum();

    // Calculate proportions
    let mut array_proportion: Vec<(Vec<u8>, f64)> = array_count
        .iter()
        .map(|(array, &count)| (array.clone(), count as f64 / total_count as f64))
        .collect();

    // Sort proportions in descending order
    array_proportion.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(Ordering::Equal));

    // Construct HapArrays
    let hap_arrays = HapArrays {
        length: snp_keys.len(),
        snp: snp_keys,
        array: array_proportion
            .into_iter()
            .map(|(array, proportion)| HapArray {
                count: array_count.get(&array).unwrap_or(&0).clone(),
                array,
                array_proportion: proportion,
                // 当前array的计数
            })
            .collect(),
        total_count,
    };

    // let elapsed_time = start_time.elapsed();
    // info!("Process_mutation_data Elapsed time: {:?}", elapsed_time);
    info!("HapArrays: {:?}", hap_arrays);
    hap_arrays
}


#[test]
pub fn test_get_read_map() {
    pretty_env_logger::init();
    let input_bam = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam";
    let snps = vec![
        ("chr5".to_string(), 70944812, 'G', 'C', 21.5, 0.4528999924659729),
        // Add more SNPs here
    ];
    let read_snp_map = get_bam_dict(input_bam, &snps);
    info!("{:?}", read_snp_map);
}





#[test]
pub fn test_get_bam_dict() {
    pretty_env_logger::init();
    let input_bam = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam";
    let snps = vec![
        ("chr5".to_string(), 70944812, 'G', 'C', 21.5, 0.4528999924659729),
        ("chr5".to_string(), 70935946, 'A', 'G', 22.649999618530273, 0.4239000082015991),
        ("chr5".to_string(), 70931073, 'A', 'C', 22.139999389648438, 0.33730000257492065),
        ("chr5".to_string(), 70940340, 'G', 'C', 21.56999969482422, 0.36230000853538513)
        // Add more SNPs here
    ];
    let read_snp_map = get_bam_dict(input_bam, &snps);
    // info!("{:?}", read_snp_map);
    let haparray = process_mutation_data(&read_snp_map);
    info!("{:?}", haparray);
}
