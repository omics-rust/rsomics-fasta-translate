use std::io::{BufWriter, Write};
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

pub fn translate_fasta(
    input: &Path,
    output: &mut dyn Write,
    frames: &[i8],
    table: u8,
) -> Result<u64> {
    if std::fs::metadata(input).is_ok_and(|m| m.len() == 0) {
        return Err(RsomicsError::InvalidInput("empty file".into()));
    }

    let mut reader = parse_fastx_file(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut count: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record.map_err(|e| RsomicsError::InvalidInput(format!("reading: {e}")))?;
        let id = std::str::from_utf8(record.id()).unwrap_or("unknown");
        let seq = record.seq();

        for &frame in frames {
            let (dna, frame_label) = if frame > 0 {
                (seq.to_vec(), format!("+{frame}"))
            } else {
                (revcomp(&seq), format!("{frame}"))
            };

            let offset = (frame.unsigned_abs() as usize).saturating_sub(1);
            if offset >= dna.len() {
                continue;
            }

            let protein = translate_seq(&dna[offset..], table);
            writeln!(out, ">{id}_frame{frame_label}").map_err(RsomicsError::Io)?;
            out.write_all(protein.as_bytes())
                .map_err(RsomicsError::Io)?;
            out.write_all(b"\n").map_err(RsomicsError::Io)?;
            count += 1;
        }
    }

    out.flush().map_err(RsomicsError::Io)?;
    Ok(count)
}

fn translate_seq(dna: &[u8], _table: u8) -> String {
    let mut protein = String::with_capacity(dna.len() / 3 + 1);
    for codon in dna.chunks(3) {
        if codon.len() < 3 {
            break;
        }
        protein.push(translate_codon(codon));
    }
    protein
}

fn translate_codon(codon: &[u8]) -> char {
    let c = [
        codon[0].to_ascii_uppercase(),
        codon[1].to_ascii_uppercase(),
        codon[2].to_ascii_uppercase(),
    ];
    match &c {
        b"TTT" | b"TTC" => 'F',
        b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
        b"ATT" | b"ATC" | b"ATA" => 'I',
        b"ATG" => 'M',
        b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
        b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
        b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
        b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
        b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
        b"TAT" | b"TAC" => 'Y',
        b"TAA" | b"TAG" | b"TGA" => '*',
        b"CAT" | b"CAC" => 'H',
        b"CAA" | b"CAG" => 'Q',
        b"AAT" | b"AAC" => 'N',
        b"AAA" | b"AAG" => 'K',
        b"GAT" | b"GAC" => 'D',
        b"GAA" | b"GAG" => 'E',
        b"TGT" | b"TGC" => 'C',
        b"TGG" => 'W',
        b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
        b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
        // Degenerate codons where all IUPAC expansions encode the same amino acid.
        // Standard IUPAC N = {A,C,G,T}; the 8 cases below are the only 4-fold
        // degenerate patterns for the standard genetic code (table 1).
        [b1, b2, b'N'] => {
            match (b1, b2) {
                (b'G', b'T') => 'V', // GTN = {GTA,GTC,GTG,GTT} all Val
                (b'G', b'C') => 'A', // GCN = {GCA,GCC,GCG,GCT} all Ala
                (b'G', b'G') => 'G', // GGN = {GGA,GGC,GGG,GGT} all Gly
                (b'C', b'T') => 'L', // CTN = {CTA,CTC,CTG,CTT} all Leu
                (b'C', b'C') => 'P', // CCN = {CCA,CCC,CCG,CCT} all Pro
                (b'C', b'G') => 'R', // CGN = {CGA,CGC,CGG,CGT} all Arg
                (b'A', b'C') => 'T', // ACN = {ACA,ACC,ACG,ACT} all Thr
                (b'T', b'C') => 'S', // TCN = {TCA,TCC,TCG,TCT} all Ser
                _ => 'X',
            }
        }
        _ => 'X',
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}
