use ndarray::prelude::*;
use ndarray::linalg::Dot;
use log::{info,debug};
use crate::bam::{get_bam_dict, process_mutation_data, HapArrays};

#[derive(Debug)]
struct Haplotype {
	array: Array1<f64>,
	hapnum: usize,
}



fn cosine_similarity(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let a_f64: Array1<f64> = a.mapv(|x| x as f64);
    let b_f64: Array1<f64> = b.mapv(|x| x as f64);
    let dot_product = a_f64.dot(&b_f64);
    let norm_a = a_f64.dot(&a_f64).sqrt();
    let norm_b = b_f64.dot(&b_f64).sqrt();
    dot_product / (norm_a * norm_b)
}



pub fn cal_haplotype(haparrays: HapArrays, snps: &Vec<(String, i64, char, char, f64, f64)>, min_proportion: f64, min_count: usize) -> usize {
	let mut hap_array = Array1::<u8>::zeros(snps.len());
	let mut hap_save = Vec::new();
	let mut array_count = 0;
	for array in haparrays.into_iter() {
		if array.array_proportion <= min_proportion || array.count < min_count {
			continue;
		}
		// debug!("{:?}", array);
		let array_add = Array1::from_vec(array.array);
		array_count += 1;
		// debug!("hap_array:{:?} array_add:{:?}", hap_array, array_add);
		hap_array = &hap_array + &array_add;
		let all_greater_than_zero = hap_array.iter().all(|&x| x > 0);
		if all_greater_than_zero {
			let hap_array_normalized = hap_array.mapv(|x| x as f64 / array_count as f64);
            hap_save.push(Haplotype { array: hap_array_normalized, hapnum: array_count });
        }
	}
    let given_vector = Array1::from_vec(snps.iter().map(|snp| snp.5 as f64).collect());
	// info!("Given vector: {:?}", given_vector);
	info!("Hap save: {:?}", hap_save);
    let mut max_similarity = -1.0;
    let mut most_similar_array = None;

    for hap in &hap_save {
        let similarity = cosine_similarity(&hap.array, &given_vector);
        if similarity > max_similarity {
            max_similarity = similarity;
            most_similar_array = Some(hap);
        }
    }

    if let Some(array) = most_similar_array {
        info!("Most similar array: {:?}; max_similarity: {:?}", array, max_similarity);
        array.hapnum
    } else {
        info!("No similar array found.");
        0
    }
	
}







// #[test]
// fn main() {
//     // 创建两个示例 Vec<f64> 向量
//     let vec1 = vec![1, 2, 3, 4];
//     let vec2 = vec![4, 5, 6, 7];

//     // 将 Vec<f64> 转换为 Array1<f64>
//     let array1 = Array1::from(vec1);
//     let array2 = Array1::from(vec2);

//     let similarity = cosine_similarity(&array1, &array2);
//     println!("Cosine similarity: {}", similarity);
// }




// #[test]
// pub fn test_hap() {
// 	pretty_env_logger::init();
// 	let input_bam = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam";
// 	let snps: Vec<(&str, i64, char, char, f64, f64)> = vec![
// 		("chr5".to_string(), 70944812, 'G', 'C', 21.5, 0.4528999924659729),
// 		("chr5".to_string(), 70935946, 'A', 'G', 22.649999618530273, 0.4239000082015991),
// 		("chr5".to_string(), 70931073, 'A', 'C', 22.139999389648438, 0.33730000257492065),
// 		("chr5".to_string(), 70940340, 'G', 'C', 21.56999969482422, 0.36230000853538513)
// 		// Add more SNPs here
// 	];
// 	let read_snp_map = get_bam_dict(input_bam, &snps);
// 	// info!("{:?}", read_snp_map);
// 	let haparrays = process_mutation_data(&read_snp_map);
// 	info!("{:?}", haparrays);
// 	cal_haplotype(haparrays, &snps);
// }
