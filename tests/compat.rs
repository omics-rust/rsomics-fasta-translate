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
