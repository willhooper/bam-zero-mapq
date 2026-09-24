use std::error::Error;
use std::ffi::{OsStr, OsString};
use std::io;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use rust_htslib::bam::{self, CompressionLevel, Read, Record};
use rust_htslib::tpool::ThreadPool;

const PROGRAM: &str = "bam-zero-unmapped-mapq";
const VERSION: &str = env!("CARGO_PKG_VERSION");

const HELP: &str = "\
bam-zero-unmapped-mapq - set MAPQ to 0 on reads carrying the BAM unmapped flag

USAGE:
    bam-zero-unmapped-mapq [OPTIONS] <INPUT.bam|-> <OUTPUT.bam|->

OPTIONS:
    -t, --threads N       BGZF worker threads [default: all logical CPUs]
    -l, --compression N   Output BAM compression level, 0-9 [default: 1]
    -q, --quiet           Do not print completion statistics
    -h, --help            Print help
    -V, --version         Print version

Use '-' for stdin or stdout. Use '--' before a path beginning with '-'.
";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output: PathBuf,
    threads: u32,
    compression: u32,
    quiet: bool,
}

enum Command {
    Run(Args),
    Help,
    Version,
}

#[derive(Debug, Default, PartialEq, Eq)]
struct Stats {
    records: u64,
    unmapped: u64,
    changed: u64,
}

fn main() -> ExitCode {
    match parse_args(std::env::args_os().skip(1)) {
        Ok(Command::Help) => {
            print!("{HELP}");
            ExitCode::SUCCESS
        }
        Ok(Command::Version) => {
            println!("{PROGRAM} {VERSION}");
            ExitCode::SUCCESS
        }
        Ok(Command::Run(args)) => {
            let quiet = args.quiet;
            let start = Instant::now();
            match run(&args) {
                Ok(stats) => {
                    if !quiet {
                        eprintln!(
                            "processed {} records; {} unmapped; {} MAPQs changed; {:.2}s",
                            stats.records,
                            stats.unmapped,
                            stats.changed,
                            start.elapsed().as_secs_f64()
                        );
                    }
                    ExitCode::SUCCESS
                }
                Err(error) => {
                    eprintln!("{PROGRAM}: error: {error}");
                    ExitCode::FAILURE
                }
            }
        }
        Err(error) => {
            eprintln!("{PROGRAM}: {error}\n\nTry '{PROGRAM} --help' for usage.");
            ExitCode::from(2)
        }
    }
}

fn parse_args<I>(args: I) -> Result<Command, String>
where
    I: IntoIterator<Item = OsString>,
{
    let mut args = args.into_iter();
    let mut positionals = Vec::with_capacity(2);
    let mut threads = default_threads();
    let mut compression = 1;
    let mut quiet = false;
    let mut options = true;

    while let Some(arg) = args.next() {
        if options && arg == OsStr::new("--") {
            options = false;
            continue;
        }

        if options {
            match arg.to_str() {
                Some("-h" | "--help") => return Ok(Command::Help),
                Some("-V" | "--version") => return Ok(Command::Version),
                Some("-q" | "--quiet") => {
                    quiet = true;
                    continue;
                }
                Some("-t" | "--threads") => {
                    let value = args
                        .next()
                        .ok_or_else(|| "--threads requires a value".to_owned())?;
                    threads = parse_threads(&value)?;
                    continue;
                }
                Some("-l" | "--compression") => {
                    let value = args
                        .next()
                        .ok_or_else(|| "--compression requires a value".to_owned())?;
                    compression = parse_compression(&value)?;
                    continue;
                }
                Some(text) if text.starts_with("--threads=") => {
                    threads = parse_threads(OsStr::new(&text[10..]))?;
                    continue;
                }
                Some(text) if text.starts_with("--compression=") => {
                    compression = parse_compression(OsStr::new(&text[14..]))?;
                    continue;
                }
                Some(text) if text.starts_with('-') && text != "-" => {
                    return Err(format!("unknown option '{text}'"));
                }
                _ => {}
            }
        }

        positionals.push(PathBuf::from(arg));
        if positionals.len() > 2 {
            return Err("expected exactly two paths: INPUT.bam and OUTPUT.bam".to_owned());
        }
    }

    if positionals.len() != 2 {
        return Err("expected exactly two paths: INPUT.bam and OUTPUT.bam".to_owned());
    }

    let output = positionals.pop().expect("length checked");
    let input = positionals.pop().expect("length checked");

    Ok(Command::Run(Args {
        input,
        output,
        threads,
        compression,
        quiet,
    }))
}

fn parse_threads(value: &OsStr) -> Result<u32, String> {
    let value = value
        .to_str()
        .ok_or_else(|| "thread count must be valid UTF-8".to_owned())?;
    let threads = value
        .parse::<u32>()
        .map_err(|_| format!("invalid thread count '{value}'"))?;
    if threads == 0 {
        return Err("thread count must be greater than zero".to_owned());
    }
    Ok(threads)
}

fn parse_compression(value: &OsStr) -> Result<u32, String> {
    let value = value
        .to_str()
        .ok_or_else(|| "compression level must be valid UTF-8".to_owned())?;
    let compression = value
        .parse::<u32>()
        .map_err(|_| format!("invalid compression level '{value}'"))?;
    if compression > 9 {
        return Err("compression level must be between 0 and 9".to_owned());
    }
    Ok(compression)
}

fn default_threads() -> u32 {
    std::thread::available_parallelism()
        .map(|count| count.get().min(u32::MAX as usize) as u32)
        .unwrap_or(1)
}

fn run(args: &Args) -> Result<Stats, Box<dyn Error>> {
    if is_file(&args.input)
        && is_file(&args.output)
        && paths_refer_to_same_file(&args.input, &args.output)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "input and output must be different files",
        )
        .into());
    }

    let mut reader = if is_stdio(&args.input) {
        bam::Reader::from_stdin().map_err(|error| contextual("cannot open stdin", error))?
    } else {
        bam::Reader::from_path(&args.input).map_err(|error| {
            contextual(
                &format!("cannot open input '{}'", args.input.display()),
                error,
            )
        })?
    };

    let header = bam::Header::from_template(reader.header());
    let mut writer = if is_stdio(&args.output) {
        bam::Writer::from_stdout(&header, bam::Format::Bam)
            .map_err(|error| contextual("cannot open stdout", error))?
    } else {
        bam::Writer::from_path(&args.output, &header, bam::Format::Bam).map_err(|error| {
            contextual(
                &format!("cannot create output '{}'", args.output.display()),
                error,
            )
        })?
    };

    writer
        .set_compression_level(CompressionLevel::Level(args.compression))
        .map_err(|error| contextual("cannot set BAM compression level", error))?;

    // Sharing one native HTSlib pool lets BGZF decompression and compression overlap
    // without oversubscribing the machine with separate input and output pools.
    let thread_pool = ThreadPool::new(args.threads)
        .map_err(|error| contextual("cannot create HTSlib thread pool", error))?;
    reader
        .set_thread_pool(&thread_pool)
        .map_err(|error| contextual("cannot attach threads to input", error))?;
    writer
        .set_thread_pool(&thread_pool)
        .map_err(|error| contextual("cannot attach threads to output", error))?;

    process_records(&mut reader, &mut writer)
}

fn process_records(
    reader: &mut bam::Reader,
    writer: &mut bam::Writer,
) -> Result<Stats, Box<dyn Error>> {
    // Reusing this allocation is significantly cheaper than `reader.records()`, which
    // creates a fresh Record for each iteration.
    let mut record = Record::new();
    let mut stats = Stats::default();

    while let Some(result) = reader.read(&mut record) {
        result.map_err(|error| {
            contextual(
                &format!("failed while reading record {}", stats.records + 1),
                error,
            )
        })?;

        stats.records += 1;
        if record.is_unmapped() {
            stats.unmapped += 1;
            if zero_unmapped_mapq(&mut record) {
                stats.changed += 1;
            }
        }

        writer.write(&record).map_err(|error| {
            contextual(
                &format!("failed while writing record {}", stats.records),
                error,
            )
        })?;
    }

    Ok(stats)
}

#[inline(always)]
fn zero_unmapped_mapq(record: &mut Record) -> bool {
    if record.mapq() == 0 {
        false
    } else {
        record.set_mapq(0);
        true
    }
}

fn contextual(context: &str, error: impl std::fmt::Display) -> io::Error {
    io::Error::new(io::ErrorKind::Other, format!("{context}: {error}"))
}

fn is_stdio(path: &Path) -> bool {
    path.as_os_str() == OsStr::new("-")
}

fn is_file(path: &Path) -> bool {
    !is_stdio(path)
}

fn paths_refer_to_same_file(left: &Path, right: &Path) -> bool {
    if left == right {
        return true;
    }

    match (std::fs::canonicalize(left), std::fs::canonicalize(right)) {
        (Ok(left), Ok(right)) => left == right,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn os_args(values: &[&str]) -> impl Iterator<Item = OsString> + '_ {
        values.iter().map(OsString::from)
    }

    #[test]
    fn unmapped_record_gets_zero_mapq() {
        let mut record = Record::new();
        record.set_unmapped();
        record.set_mapq(60);

        assert!(record.is_unmapped());
        assert!(zero_unmapped_mapq(&mut record));
        assert_eq!(record.mapq(), 0);
    }

    #[test]
    fn mapped_record_is_not_sent_to_mapq_update() {
        let mut record = Record::new();
        record.set_mapq(60);

        if record.is_unmapped() {
            zero_unmapped_mapq(&mut record);
        }

        assert_eq!(record.mapq(), 60);
    }

    #[test]
    fn unmapped_zero_mapq_is_not_counted_as_changed() {
        let mut record = Record::new();
        record.set_unmapped();

        assert!(!zero_unmapped_mapq(&mut record));
        assert_eq!(record.mapq(), 0);
    }

    #[test]
    fn parses_options_and_paths() {
        let command = parse_args(os_args(&[
            "--threads=8",
            "--compression",
            "0",
            "--quiet",
            "in.bam",
            "out.bam",
        ]))
        .expect("valid arguments");

        let Command::Run(args) = command else {
            panic!("expected run command");
        };
        assert_eq!(args.input, PathBuf::from("in.bam"));
        assert_eq!(args.output, PathBuf::from("out.bam"));
        assert_eq!(args.threads, 8);
        assert_eq!(args.compression, 0);
        assert!(args.quiet);
    }

    #[test]
    fn rejects_zero_threads() {
        let error = match parse_args(os_args(&["-t", "0", "in.bam", "out.bam"])) {
            Ok(_) => panic!("zero threads should fail"),
            Err(error) => error,
        };
        assert!(error.contains("greater than zero"));
    }

    #[test]
    fn rejects_invalid_compression() {
        let error = match parse_args(os_args(&["-l", "10", "in.bam", "out.bam"])) {
            Ok(_) => panic!("compression level 10 should fail"),
            Err(error) => error,
        };
        assert!(error.contains("between 0 and 9"));
    }
}
