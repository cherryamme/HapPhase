use clap::Parser;
use clap::builder::styling::{AnsiColor, Effects, Styles};

fn styles() -> Styles {
    Styles::styled()
        .header(AnsiColor::Yellow.on_default() | Effects::BOLD)
        .usage(AnsiColor::Yellow.on_default() | Effects::BOLD)
        .literal(AnsiColor::Blue.on_default() | Effects::BOLD)
        .placeholder(AnsiColor::Green.on_default())
}
#[command(version, author, about, long_about = None, styles = styles())]
#[command(
    help_template = "{usage-heading} {usage} \nVersion: {version} {about-section}Author:{author} Email:jiancghen2@genomics.cn/cherryamme@qq.com\n {all-args} {tab}"
)]
#[derive(Parser, Debug, Clone)]
pub struct Args {
    /// The path of input bam file
    // #[arg(short, long, default_value = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/MCGD085-12.bam")]
    #[arg(short, long)]
    pub input: String,
    /// The name of outfile
    // #[arg(short, long, default_value = "outfile.fq.gz")]
    // pub outfile: String,
    /// The path of the vcf file
    // #[arg(short, long, default_value = "/home/gushanshan/project/2024-9-30-SMN-FL-analysis/MCGD085-12/merge_output.vcf.gz")]
    #[arg(short, long)]
    pub vcf: String,
	// /// SNP regions
	// #[arg(short, long, required = false)]
	// pub snp: String,
    /// Use region
    #[arg(short, long, required = false)]
    pub region: Option<String>,
    /// het snp_num used to phase haplotype
    #[arg(short, long, default_value = "4")]
    pub snp_num: usize,
    /// min SNVratio in het snp 
    #[arg(short='m', long, default_value = "0.2")]
    pub min_ratio: f32,
    /// max SNVratio in het snp
    #[arg(short='M', long, default_value = "0.8")]
    pub max_ratio: f32,
    /// min HapArray array_proportion:
    #[arg(long, default_value = "0.1")]
    pub min_proportion: f64, 
    /// min HapArray count
    #[arg(long, default_value = "2")]
    pub min_count: usize,
    /// similarity method
    #[arg(long, default_value = "cosine",value_parser = ["cosine","cos","euclidean","euc","mse"])]
    pub method: String,

}
