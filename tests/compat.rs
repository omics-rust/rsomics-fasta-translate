use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fasta-translate"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small.fa")
}

fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

/// Per-record translated sequence (headers + line wrapping dropped).
fn seqs(fasta: &[u8]) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut started = false;
    for line in String::from_utf8_lossy(fasta).lines() {
        if line.starts_with('>') {
            if started {
                out.push(std::mem::take(&mut cur));
            }
            started = true;
        } else {
            cur.push_str(line);
        }
    }
    if started {
        out.push(cur);
    }
    out
}

// Field-level compat: the translated protein sequences must match
// `seqkit translate` (default frame 1). Ours adds a _frame+N header suffix and
// wraps differently, so compare unwrapped sequences per record, not bytes.
#[test]
fn translation_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ours = bin().arg(fixture()).output().unwrap();
    assert!(
        ours.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    let theirs = Command::new("seqkit")
        .arg("translate")
        .arg(fixture())
        .output()
        .unwrap();
    assert!(theirs.status.success());
    assert_eq!(seqs(&ours.stdout), seqs(&theirs.stdout));
}

// Degenerate codons where all IUPAC expansions encode the same amino acid must
// resolve to that amino acid, not X. seqkit translate has the same behaviour.
#[test]
fn degenerate_n_codons_match_seqkit() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    // Each of the 8 four-fold degenerate codon families that decode
    // deterministically under the standard code (table 1).
    let input = b">test\nGTNGCNGGNCTNCCNCGNACNTCN\n";
    // expected: V  A  G  L  P  R  T  S

    let tmp = std::env::temp_dir().join("rsomics-translate-degen-compat");
    std::fs::create_dir_all(&tmp).unwrap();
    let fa = tmp.join("degen.fa");
    std::fs::write(&fa, input).unwrap();

    let ours = bin().arg(&fa).output().unwrap();
    assert!(
        ours.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    let theirs = Command::new("seqkit")
        .arg("translate")
        .arg(&fa)
        .output()
        .unwrap();
    assert!(theirs.status.success());
    assert_eq!(
        seqs(&ours.stdout),
        seqs(&theirs.stdout),
        "degenerate codon translation mismatch vs seqkit"
    );
    // The result must not contain X since all codons are four-fold degenerate.
    let result = seqs(&ours.stdout);
    assert_eq!(result.len(), 1);
    assert!(
        !result[0].contains('X'),
        "expected no X in deterministic degenerate codons, got: {}",
        result[0]
    );
}

// Goldens below were captured once from `seqkit translate` (v2.9.0) so CI diffs
// ours-vs-golden per-record even where seqkit is absent.

fn golden(name: &str) -> Vec<u8> {
    std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/golden")
            .join(name),
    )
    .unwrap()
}

#[test]
fn translation_matches_seqkit_golden() {
    let ours = bin().arg(fixture()).output().unwrap();
    assert!(
        ours.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    assert_eq!(seqs(&ours.stdout), seqs(&golden("small.seqkit.faa")));
}

#[test]
fn degenerate_n_codons_match_golden() {
    let input = b">test\nGTNGCNGGNCTNCCNCGNACNTCN\n";
    let tmp = std::env::temp_dir().join("rsomics-translate-degen-golden");
    std::fs::create_dir_all(&tmp).unwrap();
    let fa = tmp.join("degen.fa");
    std::fs::write(&fa, input).unwrap();

    let ours = bin().arg(&fa).output().unwrap();
    assert!(
        ours.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    let result = seqs(&ours.stdout);
    assert_eq!(result, seqs(&golden("degen.seqkit.faa")));
    assert!(!result[0].contains('X'), "got: {}", result[0]);
}
