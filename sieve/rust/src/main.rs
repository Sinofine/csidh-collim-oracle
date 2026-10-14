use collimation_sieve_rs::rng::SplitMix64;
use collimation_sieve_rs::sieve::{
    model_queries, regularization, run, Config, LeafRounding, MapMode, ScheduleModel,
};
use collimation_sieve_rs::uint::U384;
use std::env;
use std::process::ExitCode;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

fn usage() -> &'static str {
    "usage: collimation-sieve-rs <log N> <log L> <log S> [threshold] [arity: 2|4]\n\
       [--build-map|--print-inv] [--schedule refined|paper]\n\
       [--leaf-rounding ceil|nearest] [--seed U64] [--keep-short-root]\n\n\
       [--progress NODES]\n\n\
log N = 0 selects the exact CSIDH-512 class-group order.\n\
The default refined schedule uses 2/3 for r=2 and 151/315 for r=4."
}

fn exact_csidh_order() -> U384 {
    ["3", "37", "1407181", "51593604295295867744293584889",
     "31599414504681995853008278745587832204909"]
        .iter()
        .map(|s| U384::from_decimal(s).unwrap())
        .fold(U384::ONE, |a, b| a.mul(b))
}

fn parse<T: std::str::FromStr>(s: &str, name: &str) -> Result<T, String> {
    s.parse().map_err(|_| format!("invalid {name}: {s}"))
}

fn main() -> ExitCode {
    match real_main() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => { eprintln!("error: {e}\n\n{}", usage()); ExitCode::FAILURE }
    }
}

fn real_main() -> Result<(), String> {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.len() < 3 { return Err("missing required parameters".into()); }
    let log_n: usize = parse(&args[0], "log N")?;
    let log_l: usize = parse(&args[1], "log L")?;
    let log_s: usize = parse(&args[2], "log S")?;
    if log_n > 381 { return Err("log N must not exceed 381".into()); }
    if log_l >= usize::BITS as usize || log_s >= 63 {
        return Err("L or S is too large for the host address space".into());
    }
    let mut threshold = 0.25;
    let mut arity = 2usize;
    let mut pos = 3;
    if pos < args.len() && !args[pos].starts_with("--") {
        threshold = parse(&args[pos], "threshold")?;
        pos += 1;
    }
    if pos < args.len() && !args[pos].starts_with("--") {
        arity = parse(&args[pos], "arity")?;
        pos += 1;
    }
    let mut map_mode = MapMode::PhaseOnly;
    let mut schedule_model = ScheduleModel::Refined;
    let mut leaf_rounding = LeafRounding::Ceil;
    let mut seed = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos() as u64;
    let mut retry_root = true;
    let mut progress_every = 0u128;
    while pos < args.len() {
        match args[pos].as_str() {
            "--build-map" => map_mode = MapMode::BuildMap,
            "--print-inv" => map_mode = MapMode::PrintInv,
            "--keep-short-root" => retry_root = false,
            "--schedule" => {
                pos += 1;
                schedule_model = match args.get(pos).map(String::as_str) {
                    Some("refined") => ScheduleModel::Refined,
                    Some("paper") => ScheduleModel::Paper,
                    _ => return Err("--schedule expects refined or paper".into()),
                };
            }
            "--leaf-rounding" => {
                pos += 1;
                leaf_rounding = match args.get(pos).map(String::as_str) {
                    Some("ceil") => LeafRounding::Ceil,
                    Some("nearest") => LeafRounding::Nearest,
                    _ => return Err("--leaf-rounding expects ceil or nearest".into()),
                };
            }
            "--seed" => {
                pos += 1;
                seed = parse(args.get(pos).ok_or("--seed needs a value")?, "seed")?;
            }
            "--progress" => {
                pos += 1;
                progress_every = parse(
                    args.get(pos).ok_or("--progress needs a node interval")?,
                    "progress interval",
                )?;
            }
            x => return Err(format!("unknown option: {x}")),
        }
        pos += 1;
    }
    if !(0.0..=1.0).contains(&threshold) { return Err("threshold must be in [0,1]".into()); }
    if map_mode != MapMode::PhaseOnly && arity != 4 {
        return Err("--build-map and --print-inv require arity 4".into());
    }

    let n = if log_n == 0 { exact_csidh_order() } else { U384::pow2(log_n) };
    let l = 1usize << log_l;
    let s = U384::pow2(log_s);
    let cfg = Config { n, typical_length: l, final_interval: s, arity, threshold,
                       schedule_model, leaf_rounding, map_mode, retry_root, progress_every };
    let mut rng = SplitMix64::new(seed);
    println!("model = {:?}; arity = {arity}; seed = {seed}", schedule_model);
    println!("[log N, log L, log S] = [{:.3}, {log_l}, {log_s}]", n.log2_f64());
    let started = Instant::now();
    let result = run(&cfg, &mut rng)?;
    let construction = started.elapsed();
    let print_started = Instant::now();
    if map_mode == MapMode::PrintInv {
        let map = result.map.as_ref().ok_or("--print-inv requires arity 4")?;
        for (j, x) in map.inv.iter().enumerate() {
            println!("inv[{j}] = ({},{},{},{})", x[0], x[1], x[2], x[3]);
        }
    }
    let printing = print_started.elapsed();
    let shift = rng.unit_f64() - 0.5;
    let reg = regularization(&result.phases, s, shift)?;
    let delta = result.stats.discarded as f64 / result.stats.nodes.max(1) as f64;
    let q_model = model_queries(&cfg, &result);
    let depth = result.schedule.len() - 1;
    let logs: Vec<String> = result.schedule.iter().map(|x| format!("{:.2}", x.log2_f64())).collect();
    println!("log interval schedule = [{}]", logs.join(", "));
    println!("final length = {}", result.phases.len());
    println!("construction time = {:.6} s", construction.as_secs_f64());
    if map_mode == MapMode::PrintInv { println!("inv printing time = {:.6} s", printing.as_secs_f64()); }
    println!("nodes = {}; discarded = {} ({:.2}%)", result.stats.nodes,
             result.stats.discarded, 100.0 * delta);
    println!("actual queries = {} (log2 {:.3})", result.stats.queries,
             (result.stats.queries as f64).log2());
    println!("model queries = {:.3} (log2 {:.3})", q_model, q_model.log2());
    println!("max phase-vector length = {} (log2 {:.3})", result.stats.max_length,
             (result.stats.max_length as f64).log2());
    println!("Pr[regular] = {:.3}%", 100.0 * reg.probability);
    println!("full-regularization expected bits = {:.3}", reg.expected_bits);
    println!("punctured layers = {}", reg.punctured_probabilities.len());
    println!("punctured probabilities = {:?}", reg.punctured_probabilities);
    println!("total closest probabilities [1,2,4,8] = {:?}", reg.total_closest);
    println!(
        "punctured expected bits (one closest outcome) = {:.3}",
        reg.punctured_expected_bits
    );
    println!("depth = {depth}");
    if let Some(map) = &result.map {
        println!("map entries = {}; pi entries = {}", map.inv.len(), map.pi.len());
    }
    Ok(())
}
