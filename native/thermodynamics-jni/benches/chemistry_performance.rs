use std::hint::black_box;
use std::time::{Duration, Instant};

use create_thermodynamics_jni::chemistry::dynamic::DynamicChemistryRegistry;
use create_thermodynamics_jni::chemistry::frowns::{parse_frowns, write_frowns};
use create_thermodynamics_jni::chemistry::mixture::Mixture;
use create_thermodynamics_jni::chemistry::molecule::{
    MolecularAtom, MolecularBond, MolecularStructure, ValenceSaturation,
};
use create_thermodynamics_jni::chemistry::simulation::react_for_tick;
use create_thermodynamics_jni::chemistry::substance::SubstanceId;
use create_thermodynamics_jni::chemistry::synthesis::{
    SynthesisPlanner, SynthesisReachabilityRequest, SynthesisRequest,
};
use create_thermodynamics_jni::chemistry::{
    destroy_registry_builder, destroy_registry_with_generated_reactions_builder, ChemistryRegistry,
};

struct BenchCase {
    name: &'static str,
    iterations: usize,
    total: Duration,
}

impl BenchCase {
    fn per_iteration(&self) -> Duration {
        self.total / self.iterations as u32
    }
}

fn main() {
    let filter = std::env::args().nth(1);
    let mut cases = Vec::new();
    if filter
        .as_deref()
        .is_some_and(|filter| filter.contains("synthesis"))
    {
        push_synthesis_benches(&mut cases);
        print_cases(&cases, filter.as_deref());
        return;
    }

    cases.push(run("base registry build", 50, || {
        black_box(destroy_registry_builder().unwrap().build().unwrap());
    }));

    cases.push(run("generated registry build", 10, || {
        black_box(
            destroy_registry_with_generated_reactions_builder()
                .unwrap()
                .build()
                .unwrap(),
        );
    }));

    cases.push(run("parse + canonicalize linear FROWNS", 10_000, || {
        let structure = parse_frowns(black_box("CC(C)CC(=O)O")).unwrap();
        black_box(write_frowns(&structure).unwrap());
    }));

    let benzene = parse_frowns("destroy:benzene:C,,,,,").unwrap();
    cases.push(run("canonicalize benzene", 10_000, || {
        black_box(write_frowns(black_box(&benzene)).unwrap());
    }));

    let large_cycle = large_carbon_cycle(48);
    cases.push(run("canonicalize 48 atom symmetric cycle", 250, || {
        black_box(write_frowns(black_box(&large_cycle)).unwrap());
    }));

    let base_registry = destroy_registry_builder().unwrap().build().unwrap();
    let generated_registry = destroy_registry_with_generated_reactions_builder()
        .unwrap()
        .build()
        .unwrap();
    let neutralization = neutralization_mixture(&base_registry);
    cases.push(run("reaction tick small mixture", 20_000, || {
        let mut mixture = neutralization.clone();
        black_box(react_for_tick(&base_registry, &mut mixture, 1).unwrap());
        black_box(mixture);
    }));

    let organic_mixture = organic_candidate_mixture(&generated_registry);
    cases.push(run("reaction tick generated registry", 10_000, || {
        let mut mixture = organic_mixture.clone();
        black_box(react_for_tick(&generated_registry, &mut mixture, 1).unwrap());
        black_box(mixture);
    }));

    for count in [10, 50, 100] {
        let mixture = broad_mixture(&generated_registry, count);
        cases.push(run(
            Box::leak(format!("reaction tick {count} substances").into_boxed_str()),
            5_000,
            || {
                let mut mixture = mixture.clone();
                black_box(react_for_tick(&generated_registry, &mut mixture, 1).unwrap());
                black_box(mixture);
            },
        ));
    }

    let all_substances_mixture = broad_mixture(&generated_registry, usize::MAX);
    let all_substances_count = all_substances_mixture.substances().count();
    cases.push(run(
        Box::leak(format!("reaction tick {all_substances_count} substances").into_boxed_str()),
        1_000,
        || {
            let mut mixture = all_substances_mixture.clone();
            black_box(react_for_tick(&generated_registry, &mut mixture, 1).unwrap());
            black_box(mixture);
        },
    ));

    cases.push(run("dynamic resolve cached known FROWNS", 10_000, || {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let id = registry.resolve_frowns("CC(=O)C").unwrap();
        black_box(id);
    }));

    let mut cached_dynamic_registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    cached_dynamic_registry.resolve_frowns("CCCCCCCC").unwrap();
    cases.push(run("dynamic resolve cached new FROWNS", 10_000, || {
        let id = cached_dynamic_registry.resolve_frowns("CCCCCCCC").unwrap();
        black_box(id);
    }));

    cases.push(run("dynamic generate alkene one step", 50, || {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let alkene = registry.resolve_frowns("CCCC=C").unwrap();
        black_box(registry.generate_reactions_for(&alkene, 1).unwrap());
    }));

    cases.push(run("dynamic generate methane bounded depth 2", 50, || {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let methane = registry.resolve_frowns("C").unwrap();
        black_box(registry.generate_reactions_for(&methane, 2).unwrap());
    }));

    push_synthesis_benches(&mut cases);

    print_cases(&cases, filter.as_deref());
}

fn push_synthesis_benches(cases: &mut Vec<BenchCase>) {
    let synthesis_registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    cases.push(run("synthesis hint one-step organic", 250, || {
        let routes = SynthesisPlanner::new()
            .with_max_steps(2)
            .with_max_routes(4)
            .plan_routes(
                &synthesis_registry,
                SynthesisRequest::for_substance("destroy:acetone_cyanohydrin")
                    .with_available_substance("destroy:acetone")
                    .with_available_substance("destroy:hydrogen_cyanide")
                    .with_available_substance("destroy:cyanide"),
            )
            .unwrap();
        black_box(routes);
    }));

    cases.push(run("synthesis hint deep bridgehead route", 25, || {
        let routes = SynthesisPlanner::new()
            .with_max_steps(6)
            .with_max_routes(4)
            .allow_reaction_prefix("radical_halogenation/")
            .allow_reaction_prefix("organomagnesium_formation/")
            .allow_reaction_prefix("organometallic_carboxylation/")
            .plan_routes(
                &synthesis_registry,
                SynthesisRequest::for_substance("destroy:cubanedicarboxylic_acid")
                    .with_available_substance("destroy:cubane")
                    .with_available_substance("destroy:chlorine")
                    .with_available_substance("destroy:carbon_dioxide")
                    .with_available_substance("destroy:water"),
            )
            .unwrap();
        black_box(routes);
    }));

    cases.push(run("synthesis reachability small report", 100, || {
        let report = SynthesisPlanner::new()
            .with_max_routes(2)
            .analyze_reachability(
                &synthesis_registry,
                SynthesisReachabilityRequest::new([
                    SubstanceId::from("destroy:acetone_cyanohydrin"),
                    SubstanceId::from("destroy:cubane"),
                    SubstanceId::from("destroy:iodomethane"),
                    SubstanceId::from("destroy:isopropanol"),
                ])
                .with_available_substance("destroy:propene")
                .with_available_substance("destroy:methanol")
                .with_available_substance("destroy:acetone")
                .with_available_substance("destroy:water")
                .with_available_substance("destroy:hydrogen_cyanide")
                .with_available_substance("destroy:cyanide")
                .with_available_substance("destroy:hydroiodic_acid")
                .with_available_substance("destroy:proton")
                .with_generation_iterations(1)
                .with_max_steps(4)
                .with_max_routes_per_target(2),
            )
            .unwrap();
        black_box(report);
    }));
}

fn run(name: &'static str, iterations: usize, mut f: impl FnMut()) -> BenchCase {
    for _ in 0..3 {
        f();
    }
    let started = Instant::now();
    for _ in 0..iterations {
        f();
    }
    BenchCase {
        name,
        iterations,
        total: started.elapsed(),
    }
}

fn print_cases(cases: &[BenchCase], filter: Option<&str>) {
    println!();
    println!(
        "{:<42} {:>12} {:>16} {:>16}",
        "case", "iterations", "total ms", "per op"
    );
    println!("{}", "-".repeat(90));
    for case in cases
        .iter()
        .filter(|case| filter.is_none_or(|filter| case.name.contains(filter)))
    {
        println!(
            "{:<42} {:>12} {:>16.3} {:>16}",
            case.name,
            case.iterations,
            case.total.as_secs_f64() * 1000.0,
            format_duration(case.per_iteration())
        );
    }
}

fn format_duration(duration: Duration) -> String {
    let nanos = duration.as_nanos();
    if nanos < 1_000 {
        format!("{nanos} ns")
    } else if nanos < 1_000_000 {
        format!("{:.3} us", nanos as f64 / 1_000.0)
    } else {
        format!("{:.3} ms", nanos as f64 / 1_000_000.0)
    }
}

fn neutralization_mixture(registry: &ChemistryRegistry) -> Mixture {
    let mut mixture = Mixture::new(298.0).unwrap();
    mixture
        .add_substance(registry, "destroy:proton", 0.5)
        .unwrap();
    mixture
        .add_substance(registry, "destroy:hydroxide", 0.5)
        .unwrap();
    mixture
        .add_substance(registry, "destroy:water", 1.0)
        .unwrap();
    mixture
}

fn organic_candidate_mixture(registry: &ChemistryRegistry) -> Mixture {
    let mut mixture = Mixture::new(350.0).unwrap();
    mixture
        .add_substance(registry, "destroy:ethanol", 0.5)
        .unwrap();
    mixture
        .add_substance(registry, "destroy:acetic_acid", 0.5)
        .unwrap();
    mixture
        .add_substance(registry, "destroy:chloroethane", 0.2)
        .unwrap();
    mixture
        .add_substance(registry, "destroy:hydroxide", 0.2)
        .unwrap();
    mixture
}

fn broad_mixture(registry: &ChemistryRegistry, count: usize) -> Mixture {
    let mut mixture = Mixture::new(350.0).unwrap();
    for substance in registry.substances().take(count) {
        mixture
            .add_substance(registry, substance.id.clone(), 0.01)
            .unwrap();
    }
    mixture
}

fn large_carbon_cycle(atom_count: usize) -> MolecularStructure {
    MolecularStructure {
        source_code: format!("bench:carbon-cycle-{atom_count}"),
        atoms: (0..atom_count)
            .map(|_| MolecularAtom {
                element: "C".to_string(),
                charge: 0.0,
                r_group_number: 0,
                valence_saturation: ValenceSaturation::UnsaturatedAllowed,
            })
            .collect(),
        bonds: (0..atom_count)
            .map(|index| MolecularBond {
                from: index,
                to: (index + 1) % atom_count,
                order: 1.0,
            })
            .collect(),
        stereochemistry: Vec::new(),
    }
}
