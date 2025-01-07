use log::{info,debug};
use rust_htslib::bcf::{record::Genotype, IndexedReader, Read};
use std::error::Error;

fn parse_bed(bed_path: &str) -> Result<Vec<(String, u64, u64)>, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    let mut intervals = Vec::new();
    for line in BufReader::new(File::open(bed_path)?).lines() {
        let l = line?;
        let fields: Vec<&str> = l.trim().split('\t').collect();
        if fields.len() >= 3 {
            let chr = fields[0].to_string();
            let start = fields[1].parse::<u64>()?;
            let end = fields[2].parse::<u64>()?;
            intervals.push((chr, start, end));
        }
    }
    Ok(intervals)
}
pub fn parse_vcf(file_path: &str, region: Option<String>, snp_num: usize, min_ratio: f32, max_ratio: f32, bed_file: Option<String>) -> Result<Vec<(String, i64, char, char, f64, &str, f64)>, Box<dyn Error>> {
    let mut bed_intervals = Vec::new();
    if let Some(bed_path) = bed_file {
        bed_intervals = parse_bed(&bed_path)?;
    }
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
    let mut homo_snps = Vec::new();
    let mut hetero_snps = Vec::new();
    // let mut snps = Vec::new();

    for record_result in reader.records() {
        let record = record_result?;
		// debug!("record: {:?}", record.alleles());
        let chrom = header.rid2name(record.rid().unwrap())?;
        let pos = record.pos() as i64 + 1; // VCF is 0-based, convert to 1-based
        let qual = record.qual() as f64;
        let chrom_str = std::str::from_utf8(chrom)?;
        if bed_intervals.iter().any(|(c, s, e)| c == chrom_str && (pos as u64) >= *s && (pos as u64) <= *e) {
            continue;
        }
        let alleles = record.alleles();
        let allele_count = record.allele_count();
        if allele_count != 2 || alleles[0].len() != 1 || alleles[1].len() != 1 || !record.has_filter("PASS".as_bytes()) || qual <= 5.0 {
            continue; // Skip indels and multi-allelic sites
        }
        let ref_allele = alleles[0][0] as char;
        let alt_allele = alleles[1][0] as char;
        let genotype = format!("{}",record.genotypes().unwrap().get(0));
        let expected_hetero = ["0/1", "0|1"];
        let expected_homo = ["1/1", "1|1"];
        // Check if it's a SNP (single nucleotide polymorphism)
        if ref_allele.is_alphabetic() && alt_allele.is_alphabetic() {
            if expected_hetero.contains(&genotype.as_str()) {
                if let Some(af) = record.format(b"AF").float().ok().and_then(|af| af.get(0).map(|af| af[0])) {
                    if af >= min_ratio && af <= max_ratio {
                        hetero_snps.push((std::str::from_utf8(chrom)?.to_owned(), pos, ref_allele, alt_allele, qual, "het", af as f64));
                    }
                }
            } else if expected_homo.contains(&genotype.as_str()) {
                if let Some(af) = record.format(b"AF").float().ok().and_then(|af| af.get(0).map(|af| af[0])) {
                    if af >= max_ratio {
                        homo_snps.push((std::str::from_utf8(chrom)?.to_owned(), pos, ref_allele, alt_allele, qual, "hom", af as f64));
                    }
                }
            }
		}
    }
    homo_snps.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap());
    hetero_snps.sort_by(|a, b| b.4.partial_cmp(&a.4).unwrap());

    let top_homo = homo_snps.into_iter().take(1);
    let top_hetero = hetero_snps.into_iter().take(snp_num - 1 );

    let snps: Vec<_> = top_homo.chain(top_hetero).collect();
    info!("loading SNPs: {:?}", snps);
    Ok(snps)
}

#[test]
pub fn test_parse_vcf() {
	pretty_env_logger::init();
    // let snps = parse_vcf("/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/merge_output.vcf.gz", Some("chr5:70919543-70950000".to_string()),4, 0.2 ,0.8,"test".to_string()).unwrap();
	// info!("snps: {:?}", snps);
    // assert_eq!(snps.len(), 5);
    // assert_eq!(snps[0], ("chr5".to_string(), 70935946, 'A', 'G', 22.649999618530273, 0.4239000082015991));
    // assert_eq!(snps[1], ("chr5".to_string(), 70931073, 'A', 'C', 22.139999389648438, 0.33730000257492065));
}
