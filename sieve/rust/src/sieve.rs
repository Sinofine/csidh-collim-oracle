use crate::collimation::{collimate2, collimate4, collimate4_mapped, FourListMap};
use crate::rng::SplitMix64;
use crate::uint::U384;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScheduleModel { Paper, Refined }

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LeafRounding { Ceil, Nearest }

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum MapMode { PhaseOnly, BuildMap, PrintInv }

#[derive(Debug, Default, Clone)]
pub struct Stats {
    pub queries: u128,
    pub max_length: usize,
    pub nodes: u128,
    pub discarded: u128,
}

pub struct Config {
    pub n: U384,
    pub typical_length: usize,
    pub final_interval: U384,
    pub arity: usize,
    pub threshold: f64,
    pub schedule_model: ScheduleModel,
    pub leaf_rounding: LeafRounding,
    pub map_mode: MapMode,
    pub retry_root: bool,
    pub progress_every: u128,
}

pub struct RunResult {
    pub phases: Vec<U384>,
    pub map: Option<FourListMap>,
    pub schedule: Vec<U384>,
    pub stats: Stats,
}

fn schedule_step(s: U384, l: usize, r: usize, model: ScheduleModel, n: U384) -> U384 {
    let (num, den): (f64, f64) = match (model, r) {
        (ScheduleModel::Paper, 2) => (1.0, 2.0),
        (ScheduleModel::Paper, 4) => (1.0, 4.0),
        (ScheduleModel::Refined, 2) => (2.0, 3.0),
        (ScheduleModel::Refined, 4) => (151.0, 315.0),
        _ => unreachable!(),
    };
    let raw_log = s.log2_f64() + (r - 1) as f64 * (l as f64).log2() + (num / den).log2();
    if raw_log >= n.log2_f64() {
        return n;
    }
    let lu = U384::from_u64(l as u64);
    match (model, r) {
        (ScheduleModel::Paper, 2) => s.mul(lu).div_small(2),
        (ScheduleModel::Paper, 4) => s.mul(lu).mul(lu).mul(lu).div_small(4),
        (ScheduleModel::Refined, 2) => s.mul(lu).mul_small(2).div_small(3),
        (ScheduleModel::Refined, 4) => s.mul(lu).mul(lu).mul(lu).mul_small(151).div_small(315),
        _ => unreachable!(),
    }
}

pub fn make_schedule(n: U384, s0: U384, l: usize, r: usize, model: ScheduleModel) -> Vec<U384> {
    assert!(s0 < n);
    let mut out = vec![s0];
    loop {
        let cur = *out.last().unwrap();
        let raw = schedule_step(cur, l, r, model, n);
        let next = if raw <= cur { cur.add(U384::ONE) } else { raw };
        if next >= n { break; }
        out.push(next);
    }
    out.push(n);
    out
}

pub fn run(config: &Config, rng: &mut SplitMix64) -> Result<RunResult, String> {
    if !matches!(config.arity, 2 | 4) { return Err("arity must be 2 or 4".into()); }
    if config.typical_length < 2 { return Err("L must be at least 2".into()); }
    if config.final_interval >= config.n { return Err("S must be smaller than N".into()); }
    let schedule = make_schedule(config.n, config.final_interval, config.typical_length,
                                 config.arity, config.schedule_model);
    let mut stats = Stats::default();
    let root = build_node(config, &schedule, 0, config.typical_length, true, rng, &mut stats)?;
    Ok(RunResult { phases: root.phases, map: root.map, schedule, stats })
}

struct Node { phases: Vec<U384>, map: Option<FourListMap> }

fn build_node(
    cfg: &Config, schedule: &[U384], level: usize, requested: usize, is_root: bool,
    rng: &mut SplitMix64, stats: &mut Stats,
) -> Result<Node, String> {
    loop {
        stats.nodes += 1;
        if cfg.progress_every != 0 && stats.nodes % cfg.progress_every == 0 {
            eprintln!(
                "progress: nodes={}, queries={}, discarded={}, max_length={}",
                stats.nodes, stats.queries, stats.discarded, stats.max_length
            );
        }
        let candidate = if level + 1 == schedule.len() {
            leaf(cfg.n, requested, cfg.leaf_rounding, rng, stats)?
        } else {
            internal(cfg, schedule, level, requested, is_root, rng, stats)?
        };
        stats.max_length = stats.max_length.max(candidate.phases.len());
        let keep = candidate.phases.len() as f64 >= cfg.threshold * requested as f64;
        if keep || (is_root && !cfg.retry_root) { return Ok(candidate); }
        stats.discarded += 1;
    }
}

fn leaf(n: U384, requested: usize, rounding: LeafRounding, rng: &mut SplitMix64,
        stats: &mut Stats) -> Result<Node, String> {
    let x = (requested as f64).log2();
    let ell = match rounding { LeafRounding::Ceil => x.ceil(), LeafRounding::Nearest => x.round() } as usize;
    if ell >= usize::BITS as usize { return Err("leaf vector is too large for this machine".into()); }
    let mut phases = vec![U384::ZERO];
    phases.reserve((1usize << ell).saturating_sub(1));
    for _ in 0..ell {
        let b = rng.below(n);
        let old = phases.len();
        for i in 0..old {
            let sum = phases[i].add(b);
            phases.push(if sum >= n { sum.sub(n) } else { sum });
        }
    }
    phases.sort_unstable();
    stats.queries += ell as u128;
    Ok(Node { phases, map: None })
}

fn target_product(cfg: &Config, requested: usize, s: U384, next: U384) -> f64 {
    let ratio = 2f64.powf(next.log2_f64() - s.log2_f64());
    let constant = match (cfg.schedule_model, cfg.arity) {
        (ScheduleModel::Paper, 2) => 2.0,
        (ScheduleModel::Paper, 4) => 4.0,
        (ScheduleModel::Refined, 2) => 1.5,
        (ScheduleModel::Refined, 4) => 315.0 / 151.0,
        _ => unreachable!(),
    };
    constant * requested as f64 * ratio
}

fn internal(
    cfg: &Config, schedule: &[U384], level: usize, requested: usize, is_root: bool,
    rng: &mut SplitMix64, stats: &mut Stats,
) -> Result<Node, String> {
    let target = target_product(cfg, requested, schedule[level], schedule[level + 1]);
    if !target.is_finite() { return Err("requested child product exceeds floating-point range".into()); }
    let mut children: Vec<Vec<U384>> = Vec::with_capacity(cfg.arity);
    let mut product = 1.0f64;
    for child_no in 0..cfg.arity {
        let remaining = cfg.arity - child_no;
        let wanted = (target / product).max(1.0).powf(1.0 / remaining as f64).ceil();
        if wanted > usize::MAX as f64 { return Err("requested child vector exceeds address space".into()); }
        let child = build_node(cfg, schedule, level + 1, wanted as usize, false, rng, stats)?;
        product *= child.phases.len() as f64;
        children.push(child.phases);
    }
    if children.iter().any(Vec::is_empty) { return Ok(Node { phases: Vec::new(), map: None }); }
    let sample_sum = children.iter().fold(U384::ZERO, |acc, v| acc.add(v[rng.index(v.len())]));
    let q = sample_sum.div(schedule[level]);
    if cfg.arity == 2 {
        let phases = collimate2(schedule[level], q, &children[0], &children[1]);
        Ok(Node { phases, map: None })
    } else {
        let refs = [&children[0][..], &children[1][..], &children[2][..], &children[3][..]];
        if is_root && cfg.map_mode != MapMode::PhaseOnly {
            let map = collimate4_mapped(schedule[level], q, refs);
            let phases = map.phases.clone();
            Ok(Node { phases, map: Some(map) })
        } else {
            Ok(Node { phases: collimate4(schedule[level], q, refs), map: None })
        }
    }
}

pub struct Regularization {
    pub probability: f64,
    pub expected_bits: f64,
    pub punctured_probabilities: Vec<f64>,
    pub closest_by_layer: Vec<[f64; 4]>,
    pub total_closest: [f64; 4],
    pub punctured_expected_bits: f64,
}

fn prob_theta(t: usize, values: &[U384], theta: f64) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for &b in values {
        let angle = std::f64::consts::TAU * b.to_f64() * theta;
        re += angle.cos();
        im += angle.sin();
    }
    (re * re + im * im) / (t as f64 * values.len() as f64)
}

fn closest_probabilities(t: usize, values: &[U384], shift: f64) -> [f64; 4] {
    let mut out = [0.0; 4];
    for (slot, z) in out.iter_mut().zip([1i32, 2, 4, 8]) {
        for i in 0..z {
            let theta = (shift + (i - z / 2) as f64) / t as f64;
            *slot += prob_theta(t, values, theta);
        }
    }
    out
}

pub fn regularization(phases: &[U384], s: U384, shift: f64) -> Result<Regularization, String> {
    let s_len = s.to_u64().ok_or("S is too large for numerical regularization")? as usize;
    let mut counts: BTreeMap<U384, usize> = BTreeMap::new();
    for &x in phases { *counts.entry(x).or_default() += 1; }
    let denom = phases.len().max(1) as f64;
    let minimum = if counts.len() == s_len { counts.values().copied().min().unwrap_or(0) } else { 0 };
    let p_regular = minimum as f64 * s_len as f64 / denom;
    let expected_regular = p_regular * (s_len as f64).log2();

    let max_count = counts.values().copied().max().unwrap_or(0);
    let cutoff = s_len / 64;
    let mut ps = Vec::new();
    let mut closest_by_layer = Vec::new();
    let mut total_closest = [0.0; 4];
    for layer in 1..=max_count {
        let values: Vec<U384> = counts.iter().filter_map(|(&x, &c)| (c >= layer).then_some(x)).collect();
        let len = values.len();
        if len <= cutoff { break; }
        let p = len as f64 / denom;
        ps.push(p);
        let closest = closest_probabilities(s_len, &values, shift);
        for i in 0..4 { total_closest[i] += p * closest[i]; }
        closest_by_layer.push(closest);
    }
    Ok(Regularization {
        probability: p_regular,
        expected_bits: expected_regular,
        punctured_probabilities: ps,
        closest_by_layer,
        total_closest,
        punctured_expected_bits: total_closest[0] * (s_len as f64).log2(),
    })
}

pub fn model_queries(cfg: &Config, result: &RunResult) -> f64 {
    let depth = result.schedule.len() - 1;
    let delta = result.stats.discarded as f64 / result.stats.nodes.max(1) as f64;
    let prev = result.schedule[result.schedule.len() - 2];
    let ratio = 2f64.powf(cfg.n.log2_f64() - prev.log2_f64());
    let c = match (cfg.schedule_model, cfg.arity) {
        (ScheduleModel::Paper, r) => r as f64,
        (ScheduleModel::Refined, 2) => 1.5,
        (ScheduleModel::Refined, 4) => 315.0 / 151.0,
        _ => unreachable!(),
    };
    let leaf = (c * cfg.typical_length as f64 * ratio).powf(1.0 / cfg.arity as f64);
    (cfg.arity as f64 / (1.0 - delta)).powi(depth as i32) * leaf.log2().ceil()
}
