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

// A codon whose three positions are all IUPAC gap characters translates to a
// single gap; a partial-gap codon has no consistent amino acid and is 'X'.
// Every other codon is resolved by expanding each IUPAC symbol to its base set
// (U treated as T for RNA input), translating every concrete ACGT combination
// under the standard code, and emitting the amino acid only when all agree.
fn translate_codon(codon: &[u8]) -> char {
    let c = [
        codon[0].to_ascii_uppercase(),
        codon[1].to_ascii_uppercase(),
        codon[2].to_ascii_uppercase(),
    ];
    if c == [b'-', b'-', b'-'] {
        return '-';
    }

    let (e0, e1, e2) = (iupac_bases(c[0]), iupac_bases(c[1]), iupac_bases(c[2]));
    if e0.is_empty() || e1.is_empty() || e2.is_empty() {
        return 'X';
    }

    let mut aa: Option<char> = None;
    for &b0 in e0 {
        for &b1 in e1 {
            for &b2 in e2 {
                let a = standard_codon([b0, b1, b2]);
                match aa {
                    None => aa = Some(a),
                    Some(prev) if prev != a => return 'X',
                    _ => {}
                }
            }
        }
    }
    aa.unwrap()
}

// IUPAC symbol → the set of concrete bases it denotes; empty for non-nucleotide
// characters. U aliases T so RNA translates like DNA.
fn iupac_bases(b: u8) -> &'static [u8] {
    match b {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' | b'U' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => b"",
    }
}

fn standard_codon(c: [u8; 3]) -> char {
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
        _ => 'X',
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            b'-' => b'-',
            _ => b'N',
        })
        .collect()
}
