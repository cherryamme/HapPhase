use log::{info,debug};
use rust_htslib::bcf::{IndexedReader, Read};
use std::error::Error;

pub fn parse_vcf(file_path: &str, region: Option<String>, snp_num: usize, min_ratio: f32, max_ratio: f32) -> Result<Vec<(String, i64, char, char, f64, f64)>, Box<dyn Error>> {
    let mut reader = IndexedReader::from_path(file_path).expect("Error opening VCF file");
	let header = reader.header().to_owned();
    if let Some(region_str) = region {
        let region_parts: Vec<&str> = region_str.split(':').collect();
        let rid = header.name2rid(region_parts[0].as_bytes())?;
        let positions: Vec<&str> = region_parts[1].split('-').collect();
        let start = positions[0].parse::<u64>()?;
        let end = positions.get(1).map(|&e| e.parse::<u64>().ok()).flatten();
        reader.fetch(rid, start, end)?;
    }

    let mut snps = Vec::new();

    for record_result in reader.records() {
        let record = record_result?;
		// debug!("record: {:?}", record.alleles());
        let chrom = header.rid2name(record.rid().unwrap())?;
        let pos = record.pos() as i64 + 1; // VCF is 0-based, convert to 1-based
        let qual = record.qual() as f64;

        let alleles = record.alleles();
        if alleles.len() != 2 {
            continue; // Skip indels and multi-allelic sites
        }

        let ref_allele = alleles[0][0] as char;
        let alt_allele = alleles[1][0] as char;
		let expected = ["./1", "1|1", "0/1", "0|1", "1|.", "1/1"];
        // Check if it's a SNP (single nucleotide polymorphism)
		if ref_allele.is_alphabetic() && alt_allele.is_alphabetic() {
			if let Some(af) = record.format(b"AF").float().ok().and_then(|af| af.get(0).map(|af| af[0])) {
				if af >= min_ratio && af <= max_ratio {
				snps.push((std::str::from_utf8(chrom)?.to_owned(), pos, ref_allele, alt_allele, qual, af as f64));
				}
			}
		}
    }

    // Sort by quality value in descending order and take the top 5
    snps.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap());
    snps.truncate(snp_num);

    Ok(snps)
}

#[test]
pub fn test_parse_vcf() {
	pretty_env_logger::init();
    let snps = parse_vcf("/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/merge_output.vcf.gz", Some("chr5:70919543-70950000".to_string()),4, 0.2 ,0.8).unwrap();
	info!("snps: {:?}", snps);
    // assert_eq!(snps.len(), 5);
    // assert_eq!(snps[0], ("chr5".to_string(), 70935946, 'A', 'G', 22.649999618530273, 0.4239000082015991));
    // assert_eq!(snps[1], ("chr5".to_string(), 70931073, 'A', 'C', 22.139999389648438, 0.33730000257492065));
}
