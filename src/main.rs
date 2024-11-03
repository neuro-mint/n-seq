use clap::Parser;
use colored::Colorize;
//use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::read;
use std::path::PathBuf;
use std::time::Instant;
use std::usize;

//This struct is for the cli arguements
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    ///Select sequence type:  coding | ncRNA
    #[arg(long = "seq_type")]
    seq_type: String,

    ///Number of entries to print
    #[arg(short = 'n')]
    count: usize,
}

//The following struct will be used to hold the parsed data
#[derive(Debug)]
struct SequenceRecord {
    ensemble_id: String,
    seq_type: String,
    sequence: String,
    sequence_length: usize,
    chromosome: String,
    start: i64,
    end: i64,
    ensemble_gene_id: String,
    gene_biotype: String,
    transcript_biotype: String,
}

impl SequenceRecord {
    fn print_record(&self) {
        print!(
            " Ensemble identifier      > {}\n",
            self.ensemble_id.on_black()
        );
        print!(
            " Ensemble gene identifier > {}\n",
            self.ensemble_gene_id.on_black()
        );

        if self.seq_type == "cds" {
            print!(
                " Sequence Type            > {}\n",
                "CODING".to_string().on_red()
            );
        } else if self.seq_type == "ncrna" {
            print!(
                " Sequence Type            > {}\n",
                "NON_CODING_RNA".to_string().on_red()
            );
        }

        print!(" Chromosome/Scaffold      > {} \n", self.chromosome.bold());
        print!(" Start-End                > {}-{} \n", self.start, self.end);
        print!(
            " Gene Biotype             > {} \n",
            self.gene_biotype.to_string().italic()
        );
        print!(
            " Transcript biotype       > {} \n",
            self.transcript_biotype.to_string().italic()
        );
        print!(" Sequence length          > {}\n", self.sequence_length);
        print!(
            " Molecular weight         > {} g/mol\n",
            self.molecular_weight()
        );
        print!(" GC content               > {}% \n", self.gc_content());
        print!(" AT content               > {}% \n", self.at_content());

        let freq = self.nucleotide_freq();
        print!(
            " Nucleotide frequency     > | A:{} | T:{} | G:{} | C:{} |\n",
            freq.get("A").unwrap(),
            freq.get("T").unwrap(),
            freq.get("G").unwrap(),
            freq.get("C").unwrap()
        );
        print!(" Sequence : 5' -> {} \n", self.sequence);
        print!(
            " Transcribed sequence : 5' -> {}\n",
            self.transcribe_to_rna()
        );
    }

    fn nucleotide_freq(&self) -> HashMap<&str, usize> {
        let mut nuc_freq = HashMap::new();
        nuc_freq.insert("A", self.sequence.matches("A").count());
        nuc_freq.insert("T", self.sequence.matches("T").count());
        nuc_freq.insert("G", self.sequence.matches("G").count());
        nuc_freq.insert("C", self.sequence.matches("C").count());
        return nuc_freq;
    }

    fn gc_content(&self) -> f32 {
        let gc_content: f32 = (self.sequence.matches("G").count()
            + self.sequence.matches("C").count()) as f32
            / self.sequence_length as f32
            * 100 as f32;
        return gc_content;
    }

    fn at_content(&self) -> f32 {
        let at_content: f32 = (self.sequence.matches("A").count()
            + self.sequence.matches("T").count()) as f32
            / self.sequence_length as f32
            * 100 as f32;
        return at_content;
    }

    fn molecular_weight(&self) -> f32 {
        let a: f32 = self.sequence.matches("A").count() as f32;
        let t: f32 = self.sequence.matches("T").count() as f32;
        let g: f32 = self.sequence.matches("G").count() as f32;
        let c: f32 = self.sequence.matches("C").count() as f32;
        let molecular_weight = (a * 335.2) + (t * 326.4) + (g * 351.2) + (c * 311.2) + 40.0;
        return molecular_weight;
    }

    fn transcribe_to_rna(&self) -> String {
        return self.sequence.replace("T", "U");
    }
    /*
        fn first_start_codon_position(&self) -> String {
            let seq = self.sequence.replace("T", "U");
            match seq.find("AUG") {
                Some(..) => return seq.find("AUG").unwrap().to_string(),
                None => return "no start codon present".to_string(),
            }
        }
    */
    fn restriction_sites(&self) {
        print!(" Restriction endonuclease sites \n");
        self.sau3a1_restriction_site();
        self.puviii_restriction_site();
        self.smal_restriction_site();
        print!("\n");
    }

    fn sau3a1_restriction_site(&self) {
        let x: Vec<(usize, &str)> = self.sequence.match_indices("GATC").collect();
        if x.len() == 0 {
            print!(" Sau3AI site not present\n");
        } else {
            for i in x {
                print!(" Sau3AI site at : {:?}\n", i.0);
            }
        }
    }
    fn puviii_restriction_site(&self) {
        let x: Vec<(usize, &str)> = self.sequence.match_indices("CAGCTG").collect();
        if x.len() == 0 {
            print!(" puvIII site not present\n");
        } else {
            for i in x {
                print!(" pubIII site at : {:?}\n", i.0);
            }
        }
    }
    fn smal_restriction_site(&self) {
        let x: Vec<(usize, &str)> = self.sequence.match_indices("CCCGGG").collect();
        if x.len() == 0 {
            print!(" smal*  site not present\n");
        } else {
            for i in x {
                print!(" smal*  site at : {:?}\n", i.0);
            }
        }
    }
}

fn main() {
    let start = Instant::now();
    let args = Args::parse();
    let sequence_type: &str = &args.seq_type;
    let file_path_coding = "./FASTA_files/cds.fa";
    let file_path_ncrna = "./FASTA_files/ncrna.fa";

    match sequence_type {
        "coding" => {
            let cds: Vec<SequenceRecord> = all_records(file_path_coding.into());
            for i in 0..args.count {
                cds[i].print_record();
                cds[i].restriction_sites();
            }
        }
        "ncRNA" => {
            let ncrna: Vec<SequenceRecord> = all_records(file_path_ncrna.into());
            for i in 0..args.count {
                ncrna[i].print_record();
                ncrna[i].restriction_sites();
            }
        }
        _ => print!("unnknown sequence type"),
    }

    /*
    rayon::scope(move |s| {
        s.spawn(|_s| {
            let cds: Vec<SequenceRecord> = all_records(file_path_coding.into());
            print!("Sequences in  cds file {}\n", cds.len());

            let mut protein_coding: Vec<SequenceRecord> = Vec::new();

            for i in cds {
                if i.gene_biotype == "protein_coding" {
                    protein_coding.push(i);
                } else {
                    continue;
                }
            }
            //to_dataframe(cds);
            print!(
                "Number of protein coding sequences = {}\n",
                protein_coding.len()
            );
            protein_coding[0].print_record();
        });

        s.spawn(|_s| {
            let ncrna: Vec<SequenceRecord> = all_records(file_path_ncrna.into());
            print!("Sequences in  ncrna file {}\n", ncrna.len());
            //to_dataframe(ncrna);
        });
    });
    */
    let time = start.elapsed();
    print!(" Total time to run the program program : {:?}\n", time);
}

fn all_records(file_path: PathBuf) -> Vec<SequenceRecord> {
    let start = Instant::now();
    let data: Vec<String> = read_from_fasta(file_path.into());
    let records: Vec<SequenceRecord> = split_and_parse(data);
    let time = start.elapsed();
    print!(" Total time to parse file: {:?}\n", time);
    return records;
}

//reads data from a fasta file
fn read_from_fasta(path: PathBuf) -> Vec<String> {
    return String::from_utf8(read(path).expect("Failed to read data from the FASTA file"))
        .expect("Failed to convert utf-8 to String")
        .par_split('\n')
        .map(|s: &str| s.to_string())
        .collect();
}

//parses fasta and returns a struct of the parsed data
fn split_and_parse(file_data: Vec<String>) -> Vec<SequenceRecord> {
    let mut data: (Vec<String>, Vec<String>) = (Vec::new(), Vec::new());
    let mut seq = String::new();

    //splits the contents of the FASTA file into sequences and headers
    let data_len = file_data.len();
    for i in 0..data_len {
        match i {
            0 => match file_data[i].starts_with(">") {
                true => data.0.push(file_data[i].to_owned()),
                false => continue,
            },
            _ => match file_data[i].starts_with(">") {
                true => {
                    data.0.push(file_data[i].to_owned());
                    data.1.push(seq.to_owned());
                    seq = String::new();
                }
                false => seq.push_str(&file_data[i]),
            },
        }
    }

    data.0.pop();
    assert_eq!(
        data.0.len(),
        data.1.len(),
        "The number of headers is not equal to the number of sequences"
    );

    let mut records = Vec::new();
    let dlen = data.0.len();
    for i in 0..dlen {
        let parts: Vec<String> = data.0[i].split(" ").map(|s: &str| s.to_string()).collect();
        let location: Vec<String> = parts[2].split(":").map(|s: &str| s.to_string()).collect();
        let seq: &str = &data.1[i];
        //assigns the data to the struct
        let record: SequenceRecord = SequenceRecord {
            ensemble_id: (parts[0][1..]).to_string(),
            seq_type: parts[1].to_string(),
            sequence: seq.to_string(),
            sequence_length: seq.len(),
            chromosome: location[2].to_string(),
            start: location[3].parse().unwrap(),
            end: location[4].parse().unwrap(),
            ensemble_gene_id: {
                parts[3]
                    .split(":")
                    .map(|s: &str| s.to_string())
                    .collect::<Vec<String>>()
            }[1]
            .to_owned(),
            gene_biotype: {
                parts[4]
                    .split(":")
                    .map(|s: &str| s.to_string())
                    .collect::<Vec<String>>()
            }[1]
            .to_owned(),
            transcript_biotype: {
                parts[5]
                    .split(":")
                    .map(|s: &str| s.to_string())
                    .collect::<Vec<String>>()
            }[1]
            .to_owned(),
        };
        records.push(record);
    }
    return records;
}

/*
//parses the FASTA file into a dataframe
fn to_dataframe(records: Vec<SequenceRecord>) -> DataFrame {
    let mut id: Vec<String> = Vec::new();
    let mut seq_type: Vec<String> = Vec::new();
    let mut sequence: Vec<String> = Vec::new();
    let mut transcribed_sequence: Vec<String> = Vec::new();
    let mut chromosome: Vec<String> = Vec::new();
    let mut start: Vec<i64> = Vec::new();
    let mut end: Vec<i64> = Vec::new();
    let mut ensemble_gene_id: Vec<String> = Vec::new();
    let mut gene_biotype: Vec<String> = Vec::new();
    let mut transcript_biotype: Vec<String> = Vec::new();
    let mut gc_content: Vec<f32> = Vec::new();
    let mut at_content: Vec<f32> = Vec::new();
    let mut a_freq: Vec<i8> = Vec::new();
    let mut t_freq: Vec<i8> = Vec::new();
    let mut g_freq: Vec<i8> = Vec::new();
    let mut c_freq: Vec<i8> = Vec::new();
    let mut seq_len: Vec<u8> = Vec::new();

    for x in records {
        let nuc_freq = x.nucleotide_freq();

        id.push(x.ensemble_id.to_owned());
        seq_type.push(x.seq_type.to_owned());
        sequence.push(x.sequence.to_owned());
        transcribed_sequence.push(x.transcribe_to_rna());
        chromosome.push(x.chromosome.to_owned());
        start.push(x.start.to_owned());
        end.push(x.end.to_owned());
        ensemble_gene_id.push(x.ensemble_gene_id.to_owned());
        transcript_biotype.push(x.transcript_biotype.to_owned());
        gene_biotype.push(x.gene_biotype.to_owned());
        gc_content.push(x.gc_content());
        at_content.push(x.at_content());
        a_freq.push(nuc_freq.get("A").unwrap().to_owned() as i8);
        t_freq.push(nuc_freq.get("T").unwrap().to_owned() as i8);
        g_freq.push(nuc_freq.get("G").unwrap().to_owned() as i8);
        c_freq.push(nuc_freq.get("C").unwrap().to_owned() as i8);
        seq_len.push(x.sequence_length.to_owned() as u8);
    }

    let df = df!(
        "Ensemble ID" => id,
        "Ensemble Gene ID" => ensemble_gene_id,
        "Sequence Type" => seq_type,
        "Gene biotype" => gene_biotype,
        "Transcript Biotype" => transcript_biotype,
        "Chromosome" => chromosome,
        "Start position" => start,
        "End position" => end,
        "Sequence" => sequence,
        "Sequence length" => seq_len,
        "GC content %" => gc_content,
        "AT content %" => at_content,
        //"A frequency" => a_freq,
        //"T frequency" => t_freq,
        //"G frquency" => g_freq,
        //"C frequency" => c_freq,
    )
    .unwrap();

    println!("{:?}", df.head(Some(5)));
    return df;
}*/
