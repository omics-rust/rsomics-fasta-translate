use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-fasta-translate"))
}

fn golden_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn fixture() -> PathBuf {
    golden_path("small.fa")
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

/// Per-record translated sequence (headers + line wrapping dropped): ours emits
/// a `_frame±N` header suffix and wraps differently, so compat is field-level.
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

fn golden(name: &str) -> Vec<u8> {
    std::fs::read(golden_path(name)).unwrap()
}

fn run_ok(args: &[&str], input: &Path) -> Vec<u8> {
    let ours = bin().args(args).arg(input).output().unwrap();
    assert!(
        ours.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&ours.stderr)
    );
    ours.stdout
}

// Field-level compat against a live `seqkit translate`, skipped when the oracle
// is absent. The committed-golden tests below enforce the same in CI.
#[test]
fn translation_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("skipping: seqkit not found");
        return;
    }
    let ours = run_ok(&[], &fixture());
    let theirs = Command::new("seqkit")
        .arg("translate")
        .arg(fixture())
        .output()
        .unwrap();
    assert!(theirs.status.success());
    assert_eq!(seqs(&ours), seqs(&theirs.stdout));
}

// Committed goldens captured from `seqkit translate` (v2.13.0); no live oracle.
#[test]
fn valid_dna_matches_golden() {
    let ours = run_ok(&[], &fixture());
    assert_eq!(seqs(&ours), seqs(&golden("small.seqkit.faa")));
}

#[test]
fn degenerate_four_fold_matches_golden() {
    let ours = run_ok(&[], &golden_path("degen.fa"));
    let result = seqs(&ours);
    assert_eq!(result, seqs(&golden("degen.seqkit.faa")));
    assert!(!result[0].contains('X'), "got: {}", result[0]);
}

// RNA input: U-containing codons translate as if T (seqkit: FLMG, ours had XXXG).
#[test]
fn rna_matches_golden() {
    let ours = run_ok(&[], &golden_path("rna.fa"));
    assert_eq!(seqs(&ours), seqs(&golden("rna.seqkit.faa")));
    assert_eq!(seqs(&ours), vec!["FLMG".to_string()]);
}

// IUPAC-ambiguous codons that map uniquely to one amino acid must resolve, not X
// (TTR->L, YTA->L, MGR->R, TTY->F; seqkit: LLRF, ours had XXXX).
#[test]
fn iupac_unique_matches_golden() {
    let ours = run_ok(&[], &golden_path("iupac.fa"));
    let result = seqs(&ours);
    assert_eq!(result, seqs(&golden("iupac.seqkit.faa")));
    assert_eq!(result, vec!["LLRF".to_string()]);
    assert!(!result[0].contains('X'));
}

// A full-gap '---' codon translates to a single gap (seqkit: M-*, ours had MX*).
#[test]
fn gap_codon_matches_golden() {
    let ours = run_ok(&[], &golden_path("gap.fa"));
    let result = seqs(&ours);
    assert_eq!(result, seqs(&golden("gap.seqkit.faa")));
    assert_eq!(result, vec!["M-*".to_string()]);
}

// Space-separated negative frame `-f -1` must parse and reverse-complement
// translate (seqkit: RSID/ARARAR; ours previously errored on the '-1' argument).
#[test]
fn negative_frame_space_separated_matches_golden() {
    let ours = run_ok(&["-f", "-1"], &golden_path("negframe.fa"));
    assert_eq!(seqs(&ours), seqs(&golden("negframe.m1.seqkit.faa")));
}

// Fail-loud: an empty input file exits non-zero with a stderr message.
#[test]
fn empty_input_fails_loud() {
    let dir = tempfile::tempdir().unwrap();
    let empty = dir.path().join("empty.fa");
    std::fs::write(&empty, b"").unwrap();
    let out = bin().arg(&empty).output().unwrap();
    assert!(!out.status.success());
    assert!(!out.stderr.is_empty(), "expected a stderr message");
}

// --json emits the FASTA translation followed by a single status envelope line
// that parses as one JSON document reporting success.
#[test]
fn json_envelope_is_single_ok_doc() {
    let out = bin().arg("--json").arg(fixture()).output().unwrap();
    assert!(out.status.success());
    let stdout = String::from_utf8(out.stdout).unwrap();
    let last = stdout.lines().next_back().expect("no output");
    let doc: serde_json::Value = serde_json::from_str(last).expect("envelope is not valid JSON");
    assert_eq!(doc["status"], "ok");
    assert_eq!(doc["tool"], "rsomics-fasta-translate");
}
