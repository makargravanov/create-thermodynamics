pub const DAVIES_MAX_IONIC_STRENGTH_MOLAL: f64 = 0.5;
const DAVIES_A_298_15: f64 = 0.509;

pub fn davies_log10_gamma(charge_number: i8, ionic_strength_molal: f64) -> Result<f64, ()> {
    if !ionic_strength_molal.is_finite()
        || !(0.0..=DAVIES_MAX_IONIC_STRENGTH_MOLAL).contains(&ionic_strength_molal)
    {
        return Err(());
    }

    if charge_number == 0 || ionic_strength_molal == 0.0 {
        return Ok(0.0);
    }

    let sqrt_ionic_strength = ionic_strength_molal.sqrt();
    let davies_term =
        sqrt_ionic_strength / (1.0 + sqrt_ionic_strength) - 0.3 * ionic_strength_molal;
    Ok(-DAVIES_A_298_15 * f64::from(charge_number).powi(2) * davies_term)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn davies_gamma_changes_with_ionic_strength() {
        let dilute = davies_log10_gamma(1, 0.001).unwrap();
        let stronger = davies_log10_gamma(1, 0.1).unwrap();

        assert!(stronger < dilute);
    }

    #[test]
    fn neutral_species_gamma_is_one() {
        assert_eq!(davies_log10_gamma(0, 0.3).unwrap(), 0.0);
    }

    #[test]
    fn davies_rejects_out_of_range_ionic_strength() {
        assert!(davies_log10_gamma(1, DAVIES_MAX_IONIC_STRENGTH_MOLAL + 0.001).is_err());
    }
}
