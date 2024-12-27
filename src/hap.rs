use std::fmt::Display;

use ndarray::prelude::*;
use ndarray::linalg::Dot;
use log::{info,debug};
use crate::bam::{get_bam_dict, process_mutation_data, HapArrays};

#[derive(Debug)]
struct Haplotype {
	array: Array1<f64>,
	hapnum: usize,
    hap_array_proportion: f64,
}
impl Display for Haplotype {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Haplotype: {}; Hapnum: {} proportion: {}", self.array, self.hapnum, self.hap_array_proportion)
    }
    
}

fn similarity(a: &Array1<f64>, b: &Array1<f64>, method: &str) -> f64 {
    match method {
        "cosine" | "cos" => cosine_similarity(a, b),
        "euclidean" | "euc" => euclidean_distance(a, b),
        "mse" => mean_squared_error(a, b),
        _ => panic!("Unsupported method. Use 'cosine' or 'euclidean'."),
    }
}

fn mean_squared_error(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    if a.len() != b.len() {
        panic!("Vectors must be of the same length.");
    }
    let total_error: f64 = a.iter().zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum();
    1.0 - (total_error / a.len()as f64) 
}
fn cosine_similarity(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let dot_product = a.dot(b);
    let norm_a = a.dot(a).sqrt();
    let norm_b = b.dot(b).sqrt();
    let cosine_similarity = dot_product / (norm_a * norm_b);

    // Calculate length difference weight
    let length_weight = (norm_a / norm_b).min(norm_b / norm_a);
    let adjusted_length_weight = length_weight.powf(0.5);
    // Apply the length exponent to the length weight
    let adjusted_cosine_similarityt = cosine_similarity.powf(a.len() as f64);
    // Apply the length weight to the cosine similarity
    debug!("cosine_similarity: {:?}, adjusted_cosine_similarityt:{:?}, length_weight: {:?}", cosine_similarity, adjusted_cosine_similarityt, adjusted_length_weight);
    debug!("cosine_similarity * length_weight: {}, adjusted_cosine_similarityt * length_weight: {}",cosine_similarity * adjusted_length_weight, adjusted_cosine_similarityt * adjusted_length_weight);
    adjusted_cosine_similarityt * adjusted_length_weight
}

fn euclidean_distance(a: &Array1<f64>, b: &Array1<f64>) -> f64 {
    let diff = a - b;
    diff.dot(&diff).sqrt()
}


pub fn cal_haplotype(haparrays: HapArrays, snps: &Vec<(String, i64, char, char, f64, f64)>, min_proportion: f64, min_count: usize, method: &str) -> String {
	let mut hap_array = Array1::<u8>::zeros(snps.len());
	let mut hap_save = Vec::new();
	let mut array_count = 0;
    let mut hap_array_proportion = 0.0;
    if snps.len() == 0 {
        info!("No SNP found. Maybe Haplotype num: 1.");
        return format!("1 ({} snp)",snps.len())
    }else if snps.len() == 1 {
        info!("Only one SNP found. No haplotype can be inferred.");
        return format!("0 ({} snp)",snps.len())
        
    }
	for array in haparrays.into_iter() {
		if array.array_proportion <= min_proportion || array.count < min_count {
			continue;

		}
        hap_array_proportion += array.array_proportion;
		// debug!("{:?}", array);
		let array_add = Array1::from_vec(array.array);
		array_count += 1;
		// debug!("hap_array:{:?} array_add:{:?}", hap_array, array_add);
		hap_array = &hap_array + &array_add;
		let all_greater_than_zero = hap_array.iter().all(|&x| x > 0);
		if all_greater_than_zero && hap_array_proportion >= 0.6 {
			let hap_array_normalized = hap_array.mapv(|x| x as f64 / array_count as f64);
            hap_save.push(Haplotype { array: hap_array_normalized, hapnum: array_count, hap_array_proportion: hap_array_proportion});
        }
	}
    let given_vector = Array1::from_vec(snps.iter().map(|snp| snp.5 as f64).collect());
	// info!("Given vector: {:?}", given_vector);
	// info!("Haplotype save: {:?}", hap_save);
    let mut max_similarity = -1.0;
    let mut most_similar_array = None;
    let mut haparray_with_simularity = Vec::new();
    for hap in &hap_save {
        let similarity = similarity(&hap.array, &given_vector, method);
        // debug!("Similarity: {:?}", similarity);
        haparray_with_simularity.push((hap.array.to_vec(), similarity));
        if similarity > max_similarity {
            max_similarity = similarity;
            most_similar_array = Some(hap);
        }
    }
    if let Some(array) = most_similar_array {
        info!("SNVratio array: {}", given_vector);
        info!("haparray_with_simularity: {:?}", haparray_with_simularity);
        info!("Most similar array: {}; max_similarity: {:?} ", array, max_similarity);
        if array.hap_array_proportion <=0.8 {
            info!("Not enough haplotype proportion.");
            return format!("{} ({:.2} proportion)", array.hapnum, array.hap_array_proportion)
        }
        format!("{}",array.hapnum)
    } else {
        info!("No similar array found.");
        format!("0 (No similar array)")
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
