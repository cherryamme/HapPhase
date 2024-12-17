#![allow(warnings)]
use log::{info,debug};
use clap::Parser;
pub mod bam;
mod args;
pub mod hap;
pub mod vcf;


fn main() {
    std::env::set_var("RUST_LOG", "debug");
    pretty_env_logger::init();
    let comands: Vec<String> = std::env::args().collect();
    info!("Run Command: {:?}", comands);
    let mut args = args::Args::parse();
    // args.input = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam".to_string();
    // args.vcf = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/merge_output.vcf.gz".to_string();
    let snps = vcf::parse_vcf(&args.vcf, args.region.clone(), args.snp_num, args.min_ratio, args.max_ratio).unwrap();
    // info!("snps: {:?}", snps);
    // debug!("args.input: {:?}", args.input);
    let read_snp_map = bam::get_bam_dict(&args.input, &snps);
    // info!("read_snp_map: {:?}", read_snp_map);
    let haparrays = bam::process_mutation_data(&read_snp_map);
    // info!("Hap save: {:?}", haparrays);
    let hap_num = hap::cal_haplotype(haparrays, &snps, args.min_proportion, args.min_count, &args.method);
    println!("Haplotype num: {:?}", hap_num);
}
