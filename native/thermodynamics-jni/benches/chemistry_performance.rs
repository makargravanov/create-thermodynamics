use std::hint::black_box;
use std::time::{Duration, Instant};

use create_thermodynamics_jni::chemistry::frowns::write_frowns;
use create_thermodynamics_jni::chemistry::molecule::{
    MolecularAtom, MolecularBond, MolecularStructure,
};

fn main() {
    run_torus_grid(8, 8, 3);
    run_torus_grid(10, 10, 1);
}

fn run_torus_grid(width: usize, height: usize, iterations: u32) {
    let graph = torus_grid(width, height);
    black_box(write_frowns(&graph).unwrap());

    let started = Instant::now();
    for _ in 0..iterations {
        black_box(write_frowns(&graph).unwrap());
    }
    let total = started.elapsed();

    println!();
    println!("case: canonicalize {width}x{height} carbon torus grid");
    println!("atoms: {}", graph.atom_count());
    println!("bonds: {}", graph.bond_count());
    println!("iterations: {iterations}");
    println!("total: {:.3} ms", total.as_secs_f64() * 1000.0);
    println!("per op: {}", format_duration(total / iterations));
}

fn torus_grid(width: usize, height: usize) -> MolecularStructure {
    let atom_count = width * height;
    let atoms = (0..atom_count)
        .map(|_| MolecularAtom {
            element: "C".to_string(),
            charge: 0.0,
            r_group_number: 0,
        })
        .collect();
    let mut bonds = Vec::new();
    for y in 0..height {
        for x in 0..width {
            let current = y * width + x;
            let right = y * width + ((x + 1) % width);
            let down = ((y + 1) % height) * width + x;
            if current < right || x + 1 == width {
                bonds.push(MolecularBond {
                    from: current,
                    to: right,
                    order: 1.0,
                });
            }
            if current < down || y + 1 == height {
                bonds.push(MolecularBond {
                    from: current,
                    to: down,
                    order: 1.0,
                });
            }
        }
    }
    MolecularStructure {
        source_code: format!("bench:carbon-torus-grid-{width}x{height}"),
        atoms,
        bonds,
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
