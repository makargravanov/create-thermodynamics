use super::generators::{
    bonded_hydrogens, carbonyl_atoms_from_site, first_bonded_hydrogen, generated_site_reaction_id,
};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::organic::space::SiteParticipant;
use crate::chemistry::reactive_site::ReactiveSiteKind;

#[derive(Clone)]
pub(crate) struct HalideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) halogen: usize,
    pub(crate) degree: usize,
}

#[derive(Clone)]
pub(crate) struct AlcoholSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) oxygen: usize,
    pub(crate) hydrogen: usize,
    pub(crate) degree: usize,
}

#[derive(Clone)]
pub(crate) struct AlkoxideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) oxygen: usize,
}

#[derive(Clone)]
pub(crate) struct CarbonylSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) oxygen: usize,
    pub(crate) is_ketone: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlphaCarbonylKind {
    Aldehyde,
    Ketone,
    Ester,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlphaAcidityClass {
    Ordinary,
    Activated,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlphaStericClass {
    Primary,
    Secondary,
    Tertiary,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlphaConjugation {
    None,
    Allylic,
    Benzylic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum YlideStability {
    Unstabilized,
    SemiStabilized,
    Stabilized,
}

#[derive(Clone)]
pub(crate) struct AlphaCarbonCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbonyl_carbon: usize,
    pub(crate) carbonyl_oxygen: usize,
    pub(crate) alpha_carbon: usize,
    pub(crate) alpha_hydrogens: Vec<usize>,
    pub(crate) carbonyl_kind: AlphaCarbonylKind,
    pub(crate) acidity: AlphaAcidityClass,
    pub(crate) steric_class: AlphaStericClass,
    pub(crate) conjugation: AlphaConjugation,
}

#[derive(Clone)]
pub(crate) struct CarboxylicAcidSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) hydroxyl_oxygen: usize,
    pub(crate) hydroxyl_hydrogen: usize,
}

#[derive(Clone)]
pub(crate) struct EsterSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) carbonyl_oxygen: usize,
    pub(crate) alkoxy_oxygen: usize,
}

#[derive(Clone)]
pub(crate) struct AcylChlorideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) chlorine: usize,
}

#[derive(Clone)]
pub(crate) struct AcidAnhydrideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon_a: usize,
    pub(crate) oxygen_a: usize,
    pub(crate) carbon_b: usize,
    pub(crate) oxygen_b: usize,
    pub(crate) bridge_oxygen: usize,
}

#[derive(Clone)]
pub(crate) struct AmideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) carbonyl_oxygen: usize,
    pub(crate) nitrogen: usize,
}

#[derive(Clone)]
pub(crate) struct AmineSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct ThiolSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) sulfur: usize,
    pub(crate) hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct PhosphineSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) phosphorus: usize,
}

#[derive(Clone)]
pub(crate) struct NucleophilicPhosphorusSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) phosphorus: usize,
    pub(crate) hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct PhosphoniumSaltSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) phosphorus: usize,
    pub(crate) alpha_carbon: usize,
    pub(crate) alpha_hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct PhosphorusYlideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) phosphorus: usize,
    pub(crate) alpha_carbon: usize,
    pub(crate) stability: YlideStability,
}

#[derive(Clone)]
pub(crate) struct PhosphonateCarbanionSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) phosphorus: usize,
    pub(crate) alpha_carbon: usize,
}

#[derive(Clone)]
pub(crate) struct SulfoneCarbanionSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) sulfur: usize,
    pub(crate) alpha_carbon: usize,
}

#[derive(Clone)]
pub(crate) struct SulfideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) sulfur: usize,
    pub(crate) substituent_a: usize,
    pub(crate) substituent_b: usize,
}

#[derive(Clone)]
pub(crate) struct SulfoxideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) sulfur: usize,
    pub(crate) oxygen: usize,
    pub(crate) substituent_a: usize,
    pub(crate) substituent_b: usize,
}

#[derive(Clone)]
pub(crate) struct NitrileSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) nitrogen: usize,
}

#[derive(Clone)]
pub(crate) struct NitroSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) oxygens: [usize; 2],
}

#[derive(Clone)]
pub(crate) struct OximeSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) nitrogen: usize,
    pub(crate) oxygen: usize,
    pub(crate) hydrogen: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum BisNucleophileClass {
    UreaLike,
    GuanidineLike,
    HydrazineLike,
    DiamineLike,
    AmidrazoneLike,
}

#[derive(Clone)]
pub(crate) struct HydrazoneCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) imine_nitrogen: usize,
    pub(crate) terminal_nitrogen: usize,
    pub(crate) terminal_hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct ArylHydrazoneCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) imine_nitrogen: usize,
    pub(crate) terminal_nitrogen: usize,
    pub(crate) aryl_attachment_atom: usize,
    pub(crate) terminal_hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct CyclizableHydrazoneAnnulationSite<'a> {
    pub(crate) aryl_hydrazone: ArylHydrazoneCenter<'a>,
    pub(crate) ortho_atom: usize,
    pub(crate) ortho_hydrogen: usize,
}

#[derive(Clone)]
pub(crate) struct BisNucleophileCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) first_nucleophile: usize,
    pub(crate) second_nucleophile: usize,
    pub(crate) bridge_atom: Option<usize>,
    pub(crate) class: BisNucleophileClass,
}

#[derive(Clone)]
pub(crate) struct DicarbonylElectrophileCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) first_carbonyl_carbon: usize,
    pub(crate) first_carbonyl_oxygen: usize,
    pub(crate) second_carbonyl_carbon: usize,
    pub(crate) second_carbonyl_oxygen: usize,
    pub(crate) bridge_atoms: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DicarbonylCondensationTopology {
    pub(crate) retained_carbonyl_carbon: usize,
    pub(crate) retained_carbonyl_oxygen: usize,
    pub(crate) imine_carbon: usize,
    pub(crate) imine_oxygen: usize,
    pub(crate) bridge_carbon: usize,
    pub(crate) bridge_hydrogens: Vec<usize>,
    pub(crate) imine_carbon_hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct ActivatedMethyleneCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) hydrogens: Vec<usize>,
    pub(crate) electron_withdrawing_carbons: [usize; 2],
}

#[derive(Clone)]
pub(crate) struct UreaLikeCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) hetero_atom: usize,
    pub(crate) nitrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct AmidinoCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) imine_nitrogen: usize,
    pub(crate) amino_nitrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct FormylationDonorCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) oxygen: usize,
    pub(crate) hydrogen: usize,
}

#[derive(Clone)]
pub(crate) struct UnsaturatedBondSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) high_degree_carbon: usize,
    pub(crate) low_degree_carbon: usize,
    pub(crate) is_alkyne: bool,
}

#[derive(Clone)]
pub(crate) struct BoraneSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) boron: usize,
}

#[derive(Clone)]
pub(crate) struct BorateEsterSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) oxygen: usize,
    pub(crate) boron: usize,
}

#[derive(Clone)]
pub(crate) struct IsocyanateSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) functional_carbon: usize,
    pub(crate) oxygen: usize,
}

#[derive(Clone)]
pub(crate) struct ArylHalideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) halogen: usize,
}

#[derive(Clone)]
pub(crate) struct ArylMigrationSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) ring_atoms: Vec<usize>,
    pub(crate) attachment_atoms: Vec<usize>,
}

// Protecting group center types
#[derive(Clone)]
pub(crate) struct SilylEtherCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) oxygen: usize,
    pub(crate) silicon: usize,
}

#[derive(Clone)]
pub(crate) struct AcetalCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) acetal_carbon: usize,
    pub(crate) oxygen_a: usize,
    pub(crate) oxygen_b: usize,
}

#[derive(Clone)]
pub(crate) struct BocCarbamateCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) carbonyl_carbon: usize,
    pub(crate) carbonyl_oxygen: usize,
    pub(crate) alkoxy_oxygen: usize,
    pub(crate) tert_butyl_carbon: usize,
}

#[derive(Clone)]
pub(crate) struct CbzCarbamateCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) carbonyl_carbon: usize,
    pub(crate) carbonyl_oxygen: usize,
    pub(crate) alkoxy_oxygen: usize,
}

impl<'a> SiteParticipant<'a> {
    pub(crate) fn require_kind(&self, expected: ReactiveSiteKind) -> ChemistryResult<()> {
        if self.site.kind == expected {
            Ok(())
        } else {
            Err(ChemistryError::InvalidReaction {
                reaction_id: generated_site_reaction_id("typed_site", self),
                reason: format!(
                    "expected {:?} reactive site, got {:?}",
                    expected, self.site.kind
                ),
            })
        }
    }

    pub(crate) fn halide_site(self) -> ChemistryResult<HalideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Halide)?;
        let carbon = self.site_atom_by_element("C", "halide carbon")?;
        let halogen = self
            .site
            .leaving_atom
            .or_else(|| {
                self.site.atoms.iter().copied().find(|atom| {
                    matches!(
                        self.structure.atoms[*atom].element.as_str(),
                        "F" | "Cl" | "Br" | "I"
                    )
                })
            })
            .ok_or_else(|| self.site_error("halide site has no supported halogen"))?;
        Ok(HalideSite {
            degree: self
                .site
                .substitution_degree
                .unwrap_or_else(|| self.structure.carbon_degree(carbon)),
            participant: self,
            carbon,
            halogen,
        })
    }

    pub(crate) fn alcohol_site(self) -> ChemistryResult<AlcoholSite<'a>> {
        self.require_kind(ReactiveSiteKind::Alcohol)?;
        let oxygen = self.site_atom_by_element("O", "alcohol oxygen")?;
        let carbon = self.bonded_site_atom(oxygen, "C", 1.0, "alcohol carbon")?;
        let hydrogen = first_bonded_hydrogen(self.structure, oxygen)
            .ok_or_else(|| self.site_error("alcohol oxygen has no explicit hydrogen"))?;
        Ok(AlcoholSite {
            degree: self
                .site
                .substitution_degree
                .unwrap_or_else(|| self.structure.carbon_degree(carbon)),
            participant: self,
            carbon,
            oxygen,
            hydrogen,
        })
    }

    pub(crate) fn alkoxide_site(self) -> ChemistryResult<AlkoxideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Alkoxide)?;
        let oxygen = self.site_atom_by_element("O", "alkoxide oxygen")?;
        self.bonded_site_atom(oxygen, "C", 1.0, "alkoxide carbon")?;
        Ok(AlkoxideSite {
            participant: self,
            oxygen,
        })
    }

    pub(crate) fn carbonyl_site(self) -> ChemistryResult<CarbonylSite<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl
        ) {
            return Err(self.site_error("site is not a carbonyl center"));
        }
        let (carbon, oxygen) = carbonyl_atoms_from_site(self.structure, &self.site, "carbonyl")?;
        let carbon_neighbors = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .filter(|(neighbor, order)| {
                *neighbor != oxygen
                    && self.structure.atoms[*neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(*order, 1.0)
            })
            .count();
        Ok(CarbonylSite {
            is_ketone: self.site.kind == ReactiveSiteKind::Ketone
                || (self.site.kind == ReactiveSiteKind::Carbonyl && carbon_neighbors >= 2),
            participant: self,
            carbon,
            oxygen,
        })
    }

    pub(crate) fn alpha_carbon_center(self) -> ChemistryResult<AlphaCarbonCenter<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::Enol | ReactiveSiteKind::Enolate
        ) {
            return Err(self.site_error("site is not an alpha carbon center"));
        }
        let (carbonyl_carbon, carbonyl_oxygen) =
            carbonyl_atoms_from_site(self.structure, &self.site, "alpha carbon")?;
        let alpha_carbon =
            self.site
                .atoms
                .iter()
                .copied()
                .find(|atom| {
                    *atom != carbonyl_carbon
                        && self.structure.atoms[*atom].element == "C"
                        && self.structure.neighbors(carbonyl_carbon).iter().any(
                            |(neighbor, order)| {
                                *neighbor == *atom
                                    && crate::chemistry::molecule::bond_order_matches(*order, 1.0)
                            },
                        )
                })
                .ok_or_else(|| self.site_error("alpha center has no alpha carbon"))?;
        let alpha_hydrogens = bonded_hydrogens(self.structure, alpha_carbon);
        if alpha_hydrogens.is_empty() {
            return Err(self.site_error("alpha carbon has no explicit hydrogen"));
        }
        let carbonyl_kind = alpha_carbonyl_kind(self.structure, carbonyl_carbon, carbonyl_oxygen);
        let acidity = if has_second_carbonyl_neighbor(self.structure, alpha_carbon, carbonyl_carbon)
        {
            AlphaAcidityClass::Activated
        } else {
            AlphaAcidityClass::Ordinary
        };
        let carbon_neighbors = self
            .structure
            .neighbors(alpha_carbon)
            .into_iter()
            .filter(|(neighbor, _)| self.structure.atoms[*neighbor].element == "C")
            .count();
        let steric_class = match carbon_neighbors {
            0 | 1 => AlphaStericClass::Primary,
            2 => AlphaStericClass::Secondary,
            _ => AlphaStericClass::Tertiary,
        };
        let conjugation = alpha_conjugation(self.structure, alpha_carbon, carbonyl_carbon);
        Ok(AlphaCarbonCenter {
            participant: self,
            carbonyl_carbon,
            carbonyl_oxygen,
            alpha_carbon,
            alpha_hydrogens,
            carbonyl_kind,
            acidity,
            steric_class,
            conjugation,
        })
    }

    pub(crate) fn carboxylic_acid_site(self) -> ChemistryResult<CarboxylicAcidSite<'a>> {
        self.require_kind(ReactiveSiteKind::CarboxylicAcid)?;
        let (carbon, carbonyl_oxygen) =
            carbonyl_atoms_from_site(self.structure, &self.site, "carboxylic acid")?;
        let hydroxyl_oxygen = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_oxygen
                    && self.structure.atoms[neighbor].element == "O"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("carboxylic acid has no hydroxyl oxygen"))?;
        let hydroxyl_hydrogen = first_bonded_hydrogen(self.structure, hydroxyl_oxygen)
            .ok_or_else(|| self.site_error("carboxylic acid has no explicit hydroxyl hydrogen"))?;
        Ok(CarboxylicAcidSite {
            participant: self,
            carbon,
            hydroxyl_oxygen,
            hydroxyl_hydrogen,
        })
    }

    pub(crate) fn ester_site(self) -> ChemistryResult<EsterSite<'a>> {
        self.require_kind(ReactiveSiteKind::Ester)?;
        let (carbon, carbonyl_oxygen) =
            carbonyl_atoms_from_site(self.structure, &self.site, "ester")?;
        let alkoxy_oxygen = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_oxygen
                    && self.structure.atoms[neighbor].element == "O"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("ester has no alkoxy oxygen"))?;
        Ok(EsterSite {
            participant: self,
            carbon,
            carbonyl_oxygen,
            alkoxy_oxygen,
        })
    }

    pub(crate) fn acyl_chloride_site(self) -> ChemistryResult<AcylChlorideSite<'a>> {
        self.require_kind(ReactiveSiteKind::AcylChloride)?;
        let (carbon, _) = carbonyl_atoms_from_site(self.structure, &self.site, "acyl chloride")?;
        let chlorine = self.bonded_site_atom(carbon, "Cl", 1.0, "acyl chloride chlorine")?;
        Ok(AcylChlorideSite {
            participant: self,
            carbon,
            chlorine,
        })
    }

    pub(crate) fn acid_anhydride_site(self) -> ChemistryResult<AcidAnhydrideSite<'a>> {
        self.require_kind(ReactiveSiteKind::AcidAnhydride)?;
        let carbons = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| self.structure.atoms[*atom].element == "C")
            .collect::<Vec<_>>();
        if carbons.len() != 2 {
            return Err(self.site_error("acid anhydride must contain two carbonyl carbons"));
        }
        let bridge_oxygen = self
            .site
            .atoms
            .iter()
            .copied()
            .find(|atom| {
                self.structure.atoms[*atom].element == "O"
                    && self
                        .structure
                        .neighbors(*atom)
                        .into_iter()
                        .filter(|(neighbor, order)| {
                            self.structure.atoms[*neighbor].element == "C"
                                && crate::chemistry::molecule::bond_order_matches(*order, 1.0)
                        })
                        .count()
                        == 2
            })
            .ok_or_else(|| self.site_error("acid anhydride has no bridge oxygen"))?;
        let oxygen_a = self.bonded_site_atom(carbons[0], "O", 2.0, "first carbonyl oxygen")?;
        let oxygen_b = self.bonded_site_atom(carbons[1], "O", 2.0, "second carbonyl oxygen")?;
        Ok(AcidAnhydrideSite {
            participant: self,
            carbon_a: carbons[0],
            oxygen_a,
            carbon_b: carbons[1],
            oxygen_b,
            bridge_oxygen,
        })
    }

    pub(crate) fn amide_site(self) -> ChemistryResult<AmideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Amide)?;
        let (carbon, oxygen) = carbonyl_atoms_from_site(self.structure, &self.site, "amide")?;
        let nitrogen = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != oxygen
                    && self.structure.atoms[neighbor].element == "N"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("amide has no nitrogen atom"))?;
        Ok(AmideSite {
            participant: self,
            carbon,
            carbonyl_oxygen: oxygen,
            nitrogen,
        })
    }

    pub(crate) fn amine_site(self) -> ChemistryResult<AmineSite<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::PrimaryAmine
                | ReactiveSiteKind::NonTertiaryAmine
                | ReactiveSiteKind::AmideNitrogen
        ) {
            return Err(self.site_error("site is not an amine center"));
        }
        let nitrogen = self.site_atom_by_element("N", "amine nitrogen")?;
        let hydrogens = bonded_hydrogens(self.structure, nitrogen);
        if hydrogens.is_empty() {
            return Err(self.site_error("amine has no explicit nitrogen hydrogen"));
        }
        Ok(AmineSite {
            participant: self,
            nitrogen,
            hydrogens,
        })
    }

    pub(crate) fn thiol_site(self) -> ChemistryResult<ThiolSite<'a>> {
        self.require_kind(ReactiveSiteKind::Thiol)?;
        let sulfur = self.site_atom_by_element("S", "thiol sulfur")?;
        let hydrogens = bonded_hydrogens(self.structure, sulfur);
        if hydrogens.is_empty() {
            return Err(self.site_error("thiol has no explicit sulfur hydrogen"));
        }
        Ok(ThiolSite {
            participant: self,
            sulfur,
            hydrogens,
        })
    }

    pub(crate) fn phosphine_site(self) -> ChemistryResult<PhosphineSite<'a>> {
        self.require_kind(ReactiveSiteKind::Phosphine)?;
        let phosphorus = self.site_atom_by_element("P", "phosphine phosphorus")?;
        let substituents = self
            .structure
            .neighbors(phosphorus)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                crate::chemistry::molecule::bond_order_matches(order, 1.0).then_some(neighbor)
            })
            .collect::<Vec<_>>();
        if substituents.len() != 3 {
            return Err(self.site_error("phosphine must be a neutral tertiary phosphine"));
        }
        Ok(PhosphineSite {
            participant: self,
            phosphorus,
        })
    }

    pub(crate) fn nucleophilic_phosphorus_site(
        self,
    ) -> ChemistryResult<NucleophilicPhosphorusSite<'a>> {
        self.require_kind(ReactiveSiteKind::NucleophilicPhosphorus)?;
        let phosphorus = self.site_atom_by_element("P", "nucleophilic phosphorus")?;
        let hydrogens = bonded_hydrogens(self.structure, phosphorus);
        if hydrogens.is_empty() {
            return Err(self.site_error("nucleophilic phosphorus has no P-H bond"));
        }
        Ok(NucleophilicPhosphorusSite {
            participant: self,
            phosphorus,
            hydrogens,
        })
    }

    pub(crate) fn phosphonium_salt_site(self) -> ChemistryResult<PhosphoniumSaltSite<'a>> {
        self.require_kind(ReactiveSiteKind::PhosphoniumSalt)?;
        let phosphorus = self.site_atom_by_element("P", "phosphonium phosphorus")?;
        let alpha_carbon = self
            .structure
            .neighbors(phosphorus)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("phosphonium salt has no alpha carbon"))?;
        let alpha_hydrogens = bonded_hydrogens(self.structure, alpha_carbon);
        if alpha_hydrogens.is_empty() {
            return Err(self.site_error("phosphonium salt has no explicit alpha hydrogen"));
        }
        let substituents = self
            .structure
            .neighbors(phosphorus)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                (neighbor != alpha_carbon
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .collect::<Vec<_>>();
        if substituents.len() != 3 {
            return Err(self.site_error("phosphonium salt must have three phosphorus substituents"));
        }
        Ok(PhosphoniumSaltSite {
            participant: self,
            phosphorus,
            alpha_carbon,
            alpha_hydrogens,
        })
    }

    pub(crate) fn phosphorus_ylide_site(self) -> ChemistryResult<PhosphorusYlideSite<'a>> {
        self.require_kind(ReactiveSiteKind::PhosphorusYlide)?;
        let phosphorus = self.site_atom_by_element("P", "phosphorus ylide phosphorus")?;
        let alpha_carbon = self
            .structure
            .neighbors(phosphorus)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("phosphorus ylide has no alpha carbon"))?;
        let substituents = self
            .structure
            .neighbors(phosphorus)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                (neighbor != alpha_carbon
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .collect::<Vec<_>>();
        if substituents.len() != 3 {
            return Err(self.site_error("phosphorus ylide must have three phosphorus substituents"));
        }
        let stability = phosphorus_ylide_stability(self.structure, alpha_carbon);
        Ok(PhosphorusYlideSite {
            participant: self,
            phosphorus,
            alpha_carbon,
            stability,
        })
    }

    pub(crate) fn phosphonate_carbanion_site(
        self,
    ) -> ChemistryResult<PhosphonateCarbanionSite<'a>> {
        self.require_kind(ReactiveSiteKind::PhosphonateCarbanion)?;
        let phosphorus = self.site_atom_by_element("P", "phosphonate phosphorus")?;
        let alpha_carbon = self
            .structure
            .neighbors(phosphorus)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "C"
                    && self.structure.atoms[neighbor].charge < -0.1
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("phosphonate has no anionic alpha carbon"))?;
        self.structure
            .neighbors(phosphorus)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "O"
                    && crate::chemistry::molecule::bond_order_matches(order, 2.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("phosphonate has no phosphoryl oxygen"))?;
        Ok(PhosphonateCarbanionSite {
            participant: self,
            phosphorus,
            alpha_carbon,
        })
    }

    pub(crate) fn sulfone_carbanion_site(self) -> ChemistryResult<SulfoneCarbanionSite<'a>> {
        self.require_kind(ReactiveSiteKind::SulfoneCarbanion)?;
        let sulfur = self.site_atom_by_element("S", "sulfone sulfur")?;
        let alpha_carbon = self
            .structure
            .neighbors(sulfur)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "C"
                    && self.structure.atoms[neighbor].charge < -0.1
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("sulfone has no anionic alpha carbon"))?;
        let oxygens = self
            .structure
            .neighbors(sulfur)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "O" && order >= 1.5).then_some(neighbor)
            })
            .collect::<Vec<_>>();
        if oxygens.len() < 2 {
            return Err(self.site_error("sulfone carbanion has fewer than two sulfone oxygens"));
        }
        Ok(SulfoneCarbanionSite {
            participant: self,
            sulfur,
            alpha_carbon,
        })
    }

    pub(crate) fn sulfide_site(self) -> ChemistryResult<SulfideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Sulfide)?;
        let sulfur = self.site_atom_by_element("S", "sulfide sulfur")?;
        let carbon_substituents: Vec<usize> = self
            .structure
            .neighbors(sulfur)
            .into_iter()
            .filter(|(neighbor, order)| {
                crate::chemistry::molecule::bond_order_matches(*order, 1.0)
                    && self.structure.atoms[*neighbor].element == "C"
            })
            .map(|(neighbor, _)| neighbor)
            .collect();
        if carbon_substituents.len() < 2 {
            return Err(self.site_error("sulfide has fewer than two carbon substituents"));
        }
        Ok(SulfideSite {
            participant: self,
            sulfur,
            substituent_a: carbon_substituents[0],
            substituent_b: carbon_substituents[1],
        })
    }

    pub(crate) fn sulfoxide_site(self) -> ChemistryResult<SulfoxideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Sulfoxide)?;
        let sulfur = self.site_atom_by_element("S", "sulfoxide sulfur")?;
        let oxygen = self
            .structure
            .neighbors(sulfur)
            .into_iter()
            .find(|(neighbor, order)| {
                self.structure.atoms[*neighbor].element == "O" && *order >= 1.5
            })
            .map(|(neighbor, _)| neighbor)
            .ok_or_else(|| self.site_error("sulfoxide has no S=O bond"))?;
        let carbon_substituents: Vec<usize> = self
            .structure
            .neighbors(sulfur)
            .into_iter()
            .filter(|(neighbor, order)| {
                crate::chemistry::molecule::bond_order_matches(*order, 1.0)
                    && self.structure.atoms[*neighbor].element == "C"
            })
            .map(|(neighbor, _)| neighbor)
            .collect();
        if carbon_substituents.len() < 2 {
            return Err(self.site_error("sulfoxide has fewer than two carbon substituents"));
        }
        Ok(SulfoxideSite {
            participant: self,
            sulfur,
            oxygen,
            substituent_a: carbon_substituents[0],
            substituent_b: carbon_substituents[1],
        })
    }

    pub(crate) fn nitrile_site(self) -> ChemistryResult<NitrileSite<'a>> {
        self.require_kind(ReactiveSiteKind::Nitrile)?;
        let carbon = self.site_atom_by_element("C", "nitrile carbon")?;
        let nitrogen = self.bonded_site_atom(carbon, "N", 3.0, "nitrile nitrogen")?;
        Ok(NitrileSite {
            participant: self,
            carbon,
            nitrogen,
        })
    }

    pub(crate) fn nitro_site(self) -> ChemistryResult<NitroSite<'a>> {
        self.require_kind(ReactiveSiteKind::Nitro)?;
        let nitrogen = self.site_atom_by_element("N", "nitro nitrogen")?;
        let oxygens = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| self.structure.atoms[*atom].element == "O")
            .collect::<Vec<_>>();
        let oxygens: [usize; 2] = oxygens
            .try_into()
            .map_err(|_| self.site_error("nitro center must contain exactly two oxygens"))?;
        Ok(NitroSite {
            participant: self,
            nitrogen,
            oxygens,
        })
    }

    pub(crate) fn oxime_site(self) -> ChemistryResult<OximeSite<'a>> {
        self.require_kind(ReactiveSiteKind::Oxime)?;
        let carbon = self.site_atom_by_element("C", "oxime carbon")?;
        let nitrogen = self.bonded_site_atom(carbon, "N", 2.0, "oxime nitrogen")?;
        let oxygen = self.bonded_site_atom(nitrogen, "O", 1.0, "oxime hydroxyl oxygen")?;
        let hydrogen = first_bonded_hydrogen(self.structure, oxygen)
            .ok_or_else(|| self.site_error("oxime oxygen has no explicit hydrogen"))?;
        Ok(OximeSite {
            participant: self,
            carbon,
            nitrogen,
            oxygen,
            hydrogen,
        })
    }

    pub(crate) fn hydrazone_center(self) -> ChemistryResult<HydrazoneCenter<'a>> {
        self.require_kind(ReactiveSiteKind::Hydrazone)?;
        let carbon = self.site_atom_by_element("C", "hydrazone carbon")?;
        let imine_nitrogen = self.bonded_site_atom(carbon, "N", 2.0, "hydrazone imine nitrogen")?;
        let terminal_nitrogen =
            self.bonded_site_atom(imine_nitrogen, "N", 1.0, "hydrazone terminal nitrogen")?;
        Ok(HydrazoneCenter {
            terminal_hydrogens: bonded_hydrogens(self.structure, terminal_nitrogen),
            participant: self,
            carbon,
            imine_nitrogen,
            terminal_nitrogen,
        })
    }

    pub(crate) fn aryl_hydrazone_center(self) -> ChemistryResult<ArylHydrazoneCenter<'a>> {
        let error = self.site_error("hydrazone has no aryl N-substituent");
        self.try_aryl_hydrazone_center()?.ok_or(error)
    }

    pub(crate) fn try_aryl_hydrazone_center(
        self,
    ) -> ChemistryResult<Option<ArylHydrazoneCenter<'a>>> {
        let center = self.hydrazone_center()?;
        let aryl_attachment_atom = center
            .participant
            .structure
            .neighbors(center.terminal_nitrogen)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != center.imine_nitrogen
                    && !center.terminal_hydrogens.contains(&neighbor)
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0)
                    && is_aromatic_atom(center.participant.structure, neighbor))
                .then_some(neighbor)
            })
            .map(Some)
            .unwrap_or(None);
        let Some(aryl_attachment_atom) = aryl_attachment_atom else {
            return Ok(None);
        };
        Ok(Some(ArylHydrazoneCenter {
            participant: center.participant,
            carbon: center.carbon,
            imine_nitrogen: center.imine_nitrogen,
            terminal_nitrogen: center.terminal_nitrogen,
            aryl_attachment_atom,
            terminal_hydrogens: center.terminal_hydrogens,
        }))
    }

    pub(crate) fn bis_nucleophile_center(self) -> ChemistryResult<BisNucleophileCenter<'a>> {
        let error = self.site_error("site is not a bis-nucleophile center");
        self.try_bis_nucleophile_center()?.ok_or(error)
    }

    pub(crate) fn try_bis_nucleophile_center(
        self,
    ) -> ChemistryResult<Option<BisNucleophileCenter<'a>>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::BisNucleophile | ReactiveSiteKind::UreaLike
        ) {
            return Ok(None);
        }
        let nitrogens = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| self.structure.atoms[*atom].element == "N")
            .collect::<Vec<_>>();
        if nitrogens.len() < 2 {
            return Ok(None);
        }
        let bridge_atom = self.site.atoms.iter().copied().find(|atom| {
            !nitrogens.contains(atom)
                && matches!(self.structure.atoms[*atom].element.as_str(), "C" | "N")
        });
        let class = bis_nucleophile_class(self.structure, bridge_atom, &nitrogens);
        let nucleophilic_nitrogens = if matches!(
            class,
            BisNucleophileClass::UreaLike | BisNucleophileClass::GuanidineLike
        ) {
            let Some(bridge_atom) = bridge_atom else {
                return Ok(None);
            };
            nitrogens
                .iter()
                .copied()
                .filter(|nitrogen| {
                    self.structure.neighbors(bridge_atom).into_iter().any(
                        |(neighbor, order)| {
                            neighbor == *nitrogen
                                && crate::chemistry::molecule::bond_order_matches(order, 1.0)
                        },
                    )
                })
                .collect::<Vec<_>>()
        } else {
            nitrogens
        };
        if nucleophilic_nitrogens.len() < 2 {
            return Ok(None);
        }
        Ok(Some(BisNucleophileCenter {
            participant: self,
            first_nucleophile: nucleophilic_nitrogens[0],
            second_nucleophile: nucleophilic_nitrogens[1],
            bridge_atom,
            class,
        }))
    }

    pub(crate) fn dicarbonyl_electrophile_center(
        self,
    ) -> ChemistryResult<DicarbonylElectrophileCenter<'a>> {
        self.require_kind(ReactiveSiteKind::DicarbonylElectrophile)?;
        let carbonyl_carbons = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| {
                self.structure.atoms[*atom].element == "C"
                    && self
                        .structure
                        .neighbors(*atom)
                        .into_iter()
                        .any(|(neighbor, order)| {
                            self.site.atoms.contains(&neighbor)
                                && self.structure.atoms[neighbor].element == "O"
                                && crate::chemistry::molecule::bond_order_matches(order, 2.0)
                        })
            })
            .collect::<Vec<_>>();
        if carbonyl_carbons.len() != 2 {
            return Err(self.site_error("dicarbonyl center must have two carbonyl carbons"));
        }
        let first_oxygen =
            self.bonded_site_atom(carbonyl_carbons[0], "O", 2.0, "first dicarbonyl oxygen")?;
        let second_oxygen =
            self.bonded_site_atom(carbonyl_carbons[1], "O", 2.0, "second dicarbonyl oxygen")?;
        let bridge_atoms = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| {
                !carbonyl_carbons.contains(atom)
                    && *atom != first_oxygen
                    && *atom != second_oxygen
                    && self.structure.atoms[*atom].element == "C"
            })
            .collect::<Vec<_>>();
        Ok(DicarbonylElectrophileCenter {
            participant: self,
            first_carbonyl_carbon: carbonyl_carbons[0],
            first_carbonyl_oxygen: first_oxygen,
            second_carbonyl_carbon: carbonyl_carbons[1],
            second_carbonyl_oxygen: second_oxygen,
            bridge_atoms,
        })
    }

    pub(crate) fn activated_methylene_center(
        self,
    ) -> ChemistryResult<ActivatedMethyleneCenter<'a>> {
        let center = self.dicarbonyl_electrophile_center()?;
        center.activated_methylene_center().ok_or_else(|| {
            center
                .participant
                .site_error("site is not an activated methylene center")
        })
    }

    pub(crate) fn urea_like_center(self) -> ChemistryResult<UreaLikeCenter<'a>> {
        self.require_kind(ReactiveSiteKind::UreaLike)?;
        let carbon = self.site_atom_by_element("C", "urea-like carbon")?;
        let hetero_atom = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (self.site.atoms.contains(&neighbor)
                    && matches!(
                        self.structure.atoms[neighbor].element.as_str(),
                        "O" | "S" | "N"
                    )
                    && crate::chemistry::molecule::bond_order_matches(order, 2.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("urea-like center has no double-bonded hetero atom"))?;
        let nitrogens = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                (self.site.atoms.contains(&neighbor)
                    && self.structure.atoms[neighbor].element == "N"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .collect::<Vec<_>>();
        if nitrogens.len() < 2 {
            return Err(self.site_error("urea-like center has fewer than two nitrogens"));
        }
        Ok(UreaLikeCenter {
            participant: self,
            carbon,
            hetero_atom,
            nitrogens,
        })
    }

    pub(crate) fn amidino_center(self) -> ChemistryResult<AmidinoCenter<'a>> {
        let center = self.urea_like_center()?;
        if center.participant.structure.atoms[center.hetero_atom].element != "N" {
            return Err(center
                .participant
                .site_error("amidino center must have a double-bonded nitrogen"));
        }
        Ok(AmidinoCenter {
            participant: center.participant,
            carbon: center.carbon,
            imine_nitrogen: center.hetero_atom,
            amino_nitrogens: center.nitrogens,
        })
    }

    pub(crate) fn formylation_donor_center(self) -> ChemistryResult<FormylationDonorCenter<'a>> {
        self.require_kind(ReactiveSiteKind::FormylationDonor)?;
        let carbon = self.site_atom_by_element("C", "formyl donor carbon")?;
        let oxygen = self.bonded_site_atom(carbon, "O", 2.0, "formyl donor oxygen")?;
        let hydrogen = first_bonded_hydrogen(self.structure, carbon)
            .ok_or_else(|| self.site_error("formyl donor carbon has no explicit hydrogen"))?;
        Ok(FormylationDonorCenter {
            participant: self,
            carbon,
            oxygen,
            hydrogen,
        })
    }

    pub(crate) fn unsaturated_bond_site(self) -> ChemistryResult<UnsaturatedBondSite<'a>> {
        let is_alkyne = match self.site.kind {
            ReactiveSiteKind::Alkene => false,
            ReactiveSiteKind::Alkyne => true,
            _ => return Err(self.site_error("site is not an unsaturated bond")),
        };
        let carbons = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| self.structure.atoms[*atom].element == "C")
            .collect::<Vec<_>>();
        if carbons.len() != 2 {
            return Err(self.site_error("unsaturated bond must contain exactly two carbons"));
        }
        let first_degree = self.structure.carbon_degree(carbons[0]).saturating_sub(1);
        let second_degree = self.structure.carbon_degree(carbons[1]).saturating_sub(1);
        let (high_degree_carbon, low_degree_carbon) = if second_degree > first_degree {
            (carbons[1], carbons[0])
        } else {
            (carbons[0], carbons[1])
        };
        Ok(UnsaturatedBondSite {
            participant: self,
            high_degree_carbon,
            low_degree_carbon,
            is_alkyne,
        })
    }

    pub(crate) fn borane_site(self) -> ChemistryResult<BoraneSite<'a>> {
        self.require_kind(ReactiveSiteKind::Borane)?;
        let carbon = self.site_atom_by_element("C", "borane carbon")?;
        let boron = self.bonded_site_atom(carbon, "B", 1.0, "borane boron")?;
        Ok(BoraneSite {
            participant: self,
            carbon,
            boron,
        })
    }

    pub(crate) fn borate_ester_site(self) -> ChemistryResult<BorateEsterSite<'a>> {
        self.require_kind(ReactiveSiteKind::BorateEster)?;
        let oxygen = self.site_atom_by_element("O", "borate ester oxygen")?;
        let boron = self.bonded_site_atom(oxygen, "B", 1.0, "borate ester boron")?;
        Ok(BorateEsterSite {
            participant: self,
            oxygen,
            boron,
        })
    }

    pub(crate) fn isocyanate_site(self) -> ChemistryResult<IsocyanateSite<'a>> {
        self.require_kind(ReactiveSiteKind::Isocyanate)?;
        let nitrogen = self.site_atom_by_element("N", "isocyanate nitrogen")?;
        let functional_carbon = self.bonded_site_atom(nitrogen, "C", 2.0, "isocyanate carbon")?;
        let oxygen = self.bonded_site_atom(functional_carbon, "O", 2.0, "isocyanate oxygen")?;
        Ok(IsocyanateSite {
            participant: self,
            nitrogen,
            functional_carbon,
            oxygen,
        })
    }

    pub(crate) fn aryl_halide_site(self) -> ChemistryResult<ArylHalideSite<'a>> {
        self.require_kind(ReactiveSiteKind::ArylHalide)?;
        let carbon = self.site_atom_by_element("C", "aryl halide carbon")?;
        let halogen = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (crate::chemistry::molecule::bond_order_matches(order, 1.0)
                    && matches!(
                        self.structure.atoms[neighbor].element.as_str(),
                        "F" | "Cl" | "Br" | "I"
                    ))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("aryl halide carbon has no halogen neighbor"))?;
        Ok(ArylHalideSite {
            participant: self,
            carbon,
            halogen,
        })
    }

    pub(crate) fn aryl_migration_site(self) -> ChemistryResult<ArylMigrationSite<'a>> {
        self.require_kind(ReactiveSiteKind::AromaticRing)?;
        let ring_atoms = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| is_aromatic_atom(self.structure, *atom))
            .collect::<Vec<_>>();
        if ring_atoms.is_empty() {
            return Err(self.site_error("aryl migration site has no aromatic atoms"));
        }
        let attachment_atoms = ring_atoms
            .iter()
            .copied()
            .filter(|atom| {
                self.structure
                    .neighbors(*atom)
                    .into_iter()
                    .any(|(neighbor, order)| {
                        !ring_atoms.contains(&neighbor)
                            && self.structure.atoms[neighbor].element != "H"
                            && crate::chemistry::molecule::bond_order_matches(order, 1.0)
                    })
            })
            .collect::<Vec<_>>();
        if attachment_atoms.is_empty() {
            return Err(self.site_error("aryl migration site has no external attachment atom"));
        }
        Ok(ArylMigrationSite {
            participant: self,
            ring_atoms,
            attachment_atoms,
        })
    }

    // Protecting group center methods
    pub(crate) fn silyl_ether_center(self) -> ChemistryResult<SilylEtherCenter<'a>> {
        self.require_kind(ReactiveSiteKind::SilylEther)?;
        let oxygen = self.site_atom_by_element("O", "silyl ether oxygen")?;
        let silicon = self.bonded_site_atom(oxygen, "Si", 1.0, "silyl ether silicon")?;
        self.structure
            .neighbors(oxygen)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != silicon
                    && self.structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("silyl ether oxygen has no carbon neighbor"))?;
        Ok(SilylEtherCenter {
            participant: self,
            oxygen,
            silicon,
        })
    }

    pub(crate) fn acetal_center(self) -> ChemistryResult<AcetalCenter<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::Acetal | ReactiveSiteKind::Ketal
        ) {
            return Err(self.site_error("site is not an acetal or ketal center"));
        }
        let acetal_carbon = self.site_atom_by_element("C", "acetal carbon")?;
        let oxygens = self
            .structure
            .neighbors(acetal_carbon)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                (self.structure.atoms[neighbor].element == "O"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .collect::<Vec<_>>();
        if oxygens.len() < 2 {
            return Err(self.site_error("acetal center must have two single-bonded oxygens"));
        }
        Ok(AcetalCenter {
            participant: self,
            acetal_carbon,
            oxygen_a: oxygens[0],
            oxygen_b: oxygens[1],
        })
    }

    pub(crate) fn boc_carbamate_center(self) -> ChemistryResult<BocCarbamateCenter<'a>> {
        self.require_kind(ReactiveSiteKind::BocCarbamate)?;
        let nitrogen = self.site_atom_by_element("N", "Boc carbamate nitrogen")?;
        let carbonyl_carbon =
            self.bonded_site_atom(nitrogen, "C", 1.0, "Boc carbamate carbonyl carbon")?;
        let carbonyl_oxygen =
            self.bonded_site_atom(carbonyl_carbon, "O", 2.0, "Boc carbonyl oxygen")?;
        let alkoxy_oxygen = self
            .structure
            .neighbors(carbonyl_carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_oxygen
                    && self.structure.atoms[neighbor].element == "O"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("Boc carbamate has no alkoxy oxygen"))?;
        let tert_butyl_carbon = self
            .structure
            .neighbors(alkoxy_oxygen)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_carbon
                    && self.structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("Boc carbamate has no tert-butyl carbon"))?;
        Ok(BocCarbamateCenter {
            participant: self,
            nitrogen,
            carbonyl_carbon,
            carbonyl_oxygen,
            alkoxy_oxygen,
            tert_butyl_carbon,
        })
    }

    pub(crate) fn cbz_carbamate_center(self) -> ChemistryResult<CbzCarbamateCenter<'a>> {
        self.require_kind(ReactiveSiteKind::CbzCarbamate)?;
        let nitrogen = self.site_atom_by_element("N", "Cbz carbamate nitrogen")?;
        let carbonyl_carbon =
            self.bonded_site_atom(nitrogen, "C", 1.0, "Cbz carbamate carbonyl carbon")?;
        let carbonyl_oxygen =
            self.bonded_site_atom(carbonyl_carbon, "O", 2.0, "Cbz carbonyl oxygen")?;
        let alkoxy_oxygen = self
            .structure
            .neighbors(carbonyl_carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_oxygen
                    && self.structure.atoms[neighbor].element == "O"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("Cbz carbamate has no alkoxy oxygen"))?;
        self.structure
            .neighbors(alkoxy_oxygen)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_carbon
                    && self.structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("Cbz carbamate has no benzyl carbon"))?;
        Ok(CbzCarbamateCenter {
            participant: self,
            nitrogen,
            carbonyl_carbon,
            carbonyl_oxygen,
            alkoxy_oxygen,
        })
    }

    pub(crate) fn site_atom_by_element(
        &self,
        element: &str,
        label: &str,
    ) -> ChemistryResult<usize> {
        self.site
            .atoms
            .iter()
            .copied()
            .find(|atom| self.structure.atoms[*atom].element == element)
            .ok_or_else(|| self.site_error(&format!("reactive site is missing {label}")))
    }

    pub(crate) fn bonded_site_atom(
        &self,
        parent: usize,
        element: &str,
        order: f64,
        label: &str,
    ) -> ChemistryResult<usize> {
        self.structure
            .neighbors(parent)
            .into_iter()
            .find_map(|(neighbor, bond_order)| {
                (self.site.atoms.contains(&neighbor)
                    && self.structure.atoms[neighbor].element == element
                    && crate::chemistry::molecule::bond_order_matches(bond_order, order))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error(&format!("reactive site is missing {label}")))
    }

    pub(crate) fn site_error(&self, reason: &str) -> ChemistryError {
        ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("typed_site", self),
            reason: reason.to_string(),
        }
    }
}

impl<'a> DicarbonylElectrophileCenter<'a> {
    pub(crate) fn activated_methylene_center(&self) -> Option<ActivatedMethyleneCenter<'a>> {
        if self.bridge_atoms.len() != 1 {
            return None;
        }
        let carbon = self.bridge_atoms[0];
        let hydrogens = bonded_hydrogens(self.participant.structure, carbon);
        (!hydrogens.is_empty()).then_some(ActivatedMethyleneCenter {
            participant: self.participant.clone(),
            carbon,
            hydrogens,
            electron_withdrawing_carbons: [self.first_carbonyl_carbon, self.second_carbonyl_carbon],
        })
    }

    pub(crate) fn condensation_topologies(&self) -> Vec<DicarbonylCondensationTopology> {
        if self.bridge_atoms.len() != 1 {
            return Vec::new();
        }
        let bridge_carbon = self.bridge_atoms[0];
        [
            (
                self.first_carbonyl_carbon,
                self.first_carbonyl_oxygen,
                self.second_carbonyl_carbon,
                self.second_carbonyl_oxygen,
            ),
            (
                self.second_carbonyl_carbon,
                self.second_carbonyl_oxygen,
                self.first_carbonyl_carbon,
                self.first_carbonyl_oxygen,
            ),
        ]
        .into_iter()
        .filter_map(
            |(retained_carbonyl_carbon, retained_carbonyl_oxygen, imine_carbon, imine_oxygen)| {
                let bridge_hydrogens = bonded_hydrogens(self.participant.structure, bridge_carbon);
                let imine_carbon_hydrogens =
                    bonded_hydrogens(self.participant.structure, imine_carbon);
                (!bridge_hydrogens.is_empty() && !imine_carbon_hydrogens.is_empty()).then_some(
                    DicarbonylCondensationTopology {
                        retained_carbonyl_carbon,
                        retained_carbonyl_oxygen,
                        imine_carbon,
                        imine_oxygen,
                        bridge_carbon,
                        bridge_hydrogens,
                        imine_carbon_hydrogens,
                    },
                )
            },
        )
        .collect()
    }
}

impl<'a> ArylHydrazoneCenter<'a> {
    pub(crate) fn annulation_sites(&self) -> Vec<CyclizableHydrazoneAnnulationSite<'a>> {
        if self.terminal_hydrogens.is_empty() {
            return Vec::new();
        }
        self.participant
            .structure
            .neighbors(self.aryl_attachment_atom)
            .into_iter()
            .filter_map(|(neighbor, order)| {
                (self.participant.structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.5)
                    && is_aromatic_atom(self.participant.structure, neighbor))
                .then_some(neighbor)
            })
            .filter_map(|ortho_atom| {
                first_bonded_hydrogen(self.participant.structure, ortho_atom).map(|ortho_hydrogen| {
                    CyclizableHydrazoneAnnulationSite {
                        aryl_hydrazone: self.clone(),
                        ortho_atom,
                        ortho_hydrogen,
                    }
                })
            })
            .collect()
    }
}

fn alpha_carbonyl_kind(
    structure: &crate::chemistry::molecule::MolecularStructure,
    carbonyl_carbon: usize,
    carbonyl_oxygen: usize,
) -> AlphaCarbonylKind {
    if structure
        .neighbors(carbonyl_carbon)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != carbonyl_oxygen
                && structure.atoms[neighbor].element == "O"
                && crate::chemistry::molecule::bond_order_matches(order, 1.0)
        })
    {
        AlphaCarbonylKind::Ester
    } else if first_bonded_hydrogen(structure, carbonyl_carbon).is_some() {
        AlphaCarbonylKind::Aldehyde
    } else {
        AlphaCarbonylKind::Ketone
    }
}

fn has_second_carbonyl_neighbor(
    structure: &crate::chemistry::molecule::MolecularStructure,
    alpha_carbon: usize,
    first_carbonyl: usize,
) -> bool {
    structure
        .neighbors(alpha_carbon)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != first_carbonyl
                && structure.atoms[neighbor].element == "C"
                && crate::chemistry::molecule::bond_order_matches(order, 1.0)
                && structure
                    .neighbors(neighbor)
                    .into_iter()
                    .any(|(other, bond_order)| {
                        structure.atoms[other].element == "O"
                            && crate::chemistry::molecule::bond_order_matches(bond_order, 2.0)
                    })
        })
}

fn alpha_conjugation(
    structure: &crate::chemistry::molecule::MolecularStructure,
    alpha_carbon: usize,
    carbonyl_carbon: usize,
) -> AlphaConjugation {
    for (neighbor, order) in structure.neighbors(alpha_carbon) {
        if neighbor == carbonyl_carbon || structure.atoms[neighbor].element != "C" {
            continue;
        }
        if crate::chemistry::molecule::bond_order_matches(order, 1.5) {
            return AlphaConjugation::Benzylic;
        }
        if structure
            .neighbors(neighbor)
            .into_iter()
            .any(|(other, bond_order)| {
                other != alpha_carbon
                    && structure.atoms[other].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(bond_order, 2.0)
            })
        {
            return AlphaConjugation::Allylic;
        }
    }
    AlphaConjugation::None
}

fn phosphorus_ylide_stability(
    structure: &crate::chemistry::molecule::MolecularStructure,
    alpha_carbon: usize,
) -> YlideStability {
    if structure
        .neighbors(alpha_carbon)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != alpha_carbon
                && structure.atoms[neighbor].element == "C"
                && crate::chemistry::molecule::bond_order_matches(order, 2.0)
        })
    {
        return YlideStability::Stabilized;
    }
    if structure
        .neighbors(alpha_carbon)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != alpha_carbon
                && structure.atoms[neighbor].element == "C"
                && (crate::chemistry::molecule::bond_order_matches(order, 1.5)
                    || structure
                        .neighbors(neighbor)
                        .into_iter()
                        .any(|(other, bond_order)| {
                            other != alpha_carbon
                                && structure.atoms[other].element == "C"
                                && crate::chemistry::molecule::bond_order_matches(bond_order, 2.0)
                        }))
        })
    {
        return YlideStability::SemiStabilized;
    }
    YlideStability::Unstabilized
}

fn bis_nucleophile_class(
    structure: &crate::chemistry::molecule::MolecularStructure,
    bridge_atom: Option<usize>,
    nitrogens: &[usize],
) -> BisNucleophileClass {
    let Some(bridge_atom) = bridge_atom else {
        return if nitrogens.len() >= 2 && directly_bonded(structure, nitrogens[0], nitrogens[1]) {
            BisNucleophileClass::HydrazineLike
        } else {
            BisNucleophileClass::DiamineLike
        };
    };
    if structure.atoms[bridge_atom].element == "C" {
        if structure
            .neighbors(bridge_atom)
            .into_iter()
            .any(|(neighbor, order)| {
                matches!(structure.atoms[neighbor].element.as_str(), "O" | "S")
                    && crate::chemistry::molecule::bond_order_matches(order, 2.0)
            })
        {
            return BisNucleophileClass::UreaLike;
        }
        if structure
            .neighbors(bridge_atom)
            .into_iter()
            .any(|(neighbor, order)| {
                structure.atoms[neighbor].element == "N"
                    && crate::chemistry::molecule::bond_order_matches(order, 2.0)
            })
        {
            return BisNucleophileClass::GuanidineLike;
        }
    }
    if nitrogens.len() >= 2 && directly_bonded(structure, nitrogens[0], nitrogens[1]) {
        BisNucleophileClass::HydrazineLike
    } else if structure.atoms[bridge_atom].element == "N" {
        BisNucleophileClass::AmidrazoneLike
    } else {
        BisNucleophileClass::DiamineLike
    }
}

fn directly_bonded(
    structure: &crate::chemistry::molecule::MolecularStructure,
    first: usize,
    second: usize,
) -> bool {
    structure
        .neighbors(first)
        .into_iter()
        .any(|(neighbor, _)| neighbor == second)
}

fn is_aromatic_atom(
    structure: &crate::chemistry::molecule::MolecularStructure,
    atom: usize,
) -> bool {
    structure
        .neighbors(atom)
        .into_iter()
        .any(|(_, order)| crate::chemistry::molecule::bond_order_matches(order, 1.5))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;
    use crate::chemistry::reactive_site::try_find_reactive_sites;
    use crate::chemistry::substance::Substance;

    fn participant_for(code: &str, kind: ReactiveSiteKind) -> SiteParticipant<'static> {
        let structure = Box::leak(Box::new(parse_frowns(code).unwrap()));
        let substance = Box::leak(Box::new({
            let mut substance =
                Substance::new("test:typed_center", 0, 1.0, 1000.0, 500.0, 100.0, 20_000.0);
            substance.molecular_structure = Some(structure.clone());
            substance
        }));
        let site = try_find_reactive_sites(structure)
            .unwrap()
            .into_iter()
            .find(|site| site.kind == kind)
            .unwrap_or_else(|| panic!("missing reactive site {kind:?}"));
        SiteParticipant {
            substance,
            structure,
            site,
        }
    }

    #[test]
    fn hydrazone_center_is_typed_from_graph_atoms() {
        let center = participant_for(
            "destroy:graph:atoms=C.N.N.H.H.H.H;\
             bonds=0-d-1,1-s-2,0-s-3,0-s-4,2-s-5,2-s-6",
            ReactiveSiteKind::Hydrazone,
        )
        .hydrazone_center()
        .unwrap();
        assert_eq!(
            center.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(center.carbon, 0);
        assert_eq!(center.imine_nitrogen, 1);
        assert_eq!(center.terminal_nitrogen, 2);
        assert_eq!(center.terminal_hydrogens.len(), 2);
    }

    #[test]
    fn urea_like_center_is_a_bis_nucleophile() {
        let participant = participant_for(
            "destroy:graph:atoms=C.O.N.N.H.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,2-s-4,2-s-5,3-s-6,3-s-7",
            ReactiveSiteKind::UreaLike,
        );
        let urea = participant.clone().urea_like_center().unwrap();
        assert_eq!(urea.participant.substance.id.as_str(), "test:typed_center");
        assert_eq!(urea.carbon, 0);
        assert_eq!(urea.hetero_atom, 1);
        assert_eq!(urea.nitrogens.len(), 2);

        let bis_nucleophile = participant.bis_nucleophile_center().unwrap();
        assert_eq!(
            bis_nucleophile.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(bis_nucleophile.first_nucleophile, 2);
        assert_eq!(bis_nucleophile.second_nucleophile, 3);
        assert_eq!(bis_nucleophile.bridge_atom, Some(0));
        assert_eq!(bis_nucleophile.class, BisNucleophileClass::UreaLike);
    }

    #[test]
    fn guanidine_like_center_uses_single_bonded_amino_nitrogens() {
        let bis_nucleophile = participant_for(
            "destroy:graph:atoms=C.N.N.N.H.H.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,1-s-4,2-s-5,2-s-6,3-s-7,3-s-8",
            ReactiveSiteKind::UreaLike,
        )
        .bis_nucleophile_center()
        .unwrap();
        assert_eq!(bis_nucleophile.class, BisNucleophileClass::GuanidineLike);
        assert_eq!(bis_nucleophile.first_nucleophile, 2);
        assert_eq!(bis_nucleophile.second_nucleophile, 3);
        assert_eq!(bis_nucleophile.bridge_atom, Some(0));
    }

    #[test]
    fn dicarbonyl_and_formyl_donor_centers_are_typed() {
        let dicarbonyl = participant_for(
            "destroy:graph:atoms=C.O.C.C.O.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-d-4,0-s-5,2-s-6,2-s-7,3-s-8",
            ReactiveSiteKind::DicarbonylElectrophile,
        )
        .dicarbonyl_electrophile_center()
        .unwrap();
        assert_eq!(
            dicarbonyl.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(dicarbonyl.first_carbonyl_carbon, 0);
        assert_eq!(dicarbonyl.first_carbonyl_oxygen, 1);
        assert_eq!(dicarbonyl.second_carbonyl_carbon, 3);
        assert_eq!(dicarbonyl.second_carbonyl_oxygen, 4);
        assert_eq!(dicarbonyl.bridge_atoms, vec![2]);

        let formyl = participant_for(
            "destroy:graph:atoms=C.O.N.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,2-s-4,2-s-5",
            ReactiveSiteKind::FormylationDonor,
        )
        .formylation_donor_center()
        .unwrap();
        assert_eq!(
            formyl.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(formyl.carbon, 0);
        assert_eq!(formyl.oxygen, 1);
        assert_eq!(formyl.hydrogen, 3);
    }

    #[test]
    fn aryl_hydrazone_and_aryl_migration_centers_are_typed() {
        let code = "destroy:graph:atoms=C.C.C.C.C.C.N.N.C.H.H.H.H.H.H.H.H;\
             bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,\
             0-s-6,6-s-7,7-d-8,\
             1-s-9,2-s-10,3-s-11,4-s-12,5-s-13,6-s-14,8-s-15,8-s-16";
        let aryl_hydrazone = participant_for(code, ReactiveSiteKind::Hydrazone)
            .aryl_hydrazone_center()
            .unwrap();
        assert_eq!(
            aryl_hydrazone.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(aryl_hydrazone.carbon, 8);
        assert_eq!(aryl_hydrazone.imine_nitrogen, 7);
        assert_eq!(aryl_hydrazone.terminal_nitrogen, 6);
        assert_eq!(aryl_hydrazone.aryl_attachment_atom, 0);
        assert_eq!(aryl_hydrazone.terminal_hydrogens, vec![14]);

        let migration = participant_for(code, ReactiveSiteKind::AromaticRing)
            .aryl_migration_site()
            .unwrap();
        assert_eq!(
            migration.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(migration.ring_atoms.len(), 6);
        assert_eq!(migration.attachment_atoms, vec![0]);
    }

    #[test]
    fn aryl_hydrazone_annulation_sites_require_explicit_ortho_hydrogens() {
        let code = "destroy:graph:atoms=C.C.C.C.C.C.N.N.C.H.H.H.H.H.H.H.H;\
             bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,\
             0-s-6,6-s-7,7-d-8,\
             1-s-9,2-s-10,3-s-11,4-s-12,5-s-13,6-s-14,8-s-15,8-s-16";
        let aryl_hydrazone = participant_for(code, ReactiveSiteKind::Hydrazone)
            .aryl_hydrazone_center()
            .unwrap();
        let sites = aryl_hydrazone.annulation_sites();
        let ortho_atoms = sites.iter().map(|site| site.ortho_atom).collect::<Vec<_>>();
        let ortho_hydrogens = sites
            .iter()
            .map(|site| site.ortho_hydrogen)
            .collect::<Vec<_>>();
        assert_eq!(ortho_atoms, vec![1, 5]);
        assert_eq!(ortho_hydrogens, vec![9, 13]);
        assert!(sites
            .iter()
            .all(|site| site.aryl_hydrazone.carbon == aryl_hydrazone.carbon));
    }

    #[test]
    fn aryl_hydrazone_annulation_skips_blocked_ortho_position() {
        let code = "destroy:graph:atoms=C.C.C.C.C.C.N.N.C.C.H.H.H.H.H.H.H.H.H;\
             bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,\
             0-s-6,6-s-7,7-d-8,1-s-9,\
             2-s-10,3-s-11,4-s-12,5-s-13,6-s-14,8-s-15,8-s-16,9-s-17,9-s-18";
        let aryl_hydrazone = participant_for(code, ReactiveSiteKind::Hydrazone)
            .aryl_hydrazone_center()
            .unwrap();
        let sites = aryl_hydrazone.annulation_sites();
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].ortho_atom, 5);
        assert_eq!(sites[0].ortho_hydrogen, 13);
    }

    #[test]
    fn activated_methylene_and_amidino_centers_are_typed() {
        let activated_methylene = participant_for(
            "destroy:graph:atoms=C.O.C.C.O.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-d-4,0-s-5,2-s-6,2-s-7,3-s-8",
            ReactiveSiteKind::DicarbonylElectrophile,
        )
        .activated_methylene_center()
        .unwrap();
        assert_eq!(
            activated_methylene.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(activated_methylene.carbon, 2);
        assert_eq!(activated_methylene.hydrogens, vec![6, 7]);
        assert_eq!(activated_methylene.electron_withdrawing_carbons, [0, 3]);

        let amidino = participant_for(
            "destroy:graph:atoms=C.N.N.N.H.H.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,1-s-4,2-s-5,2-s-6,3-s-7,3-s-8",
            ReactiveSiteKind::UreaLike,
        )
        .amidino_center()
        .unwrap();
        assert_eq!(
            amidino.participant.substance.id.as_str(),
            "test:typed_center"
        );
        assert_eq!(amidino.carbon, 0);
        assert_eq!(amidino.imine_nitrogen, 1);
        assert_eq!(amidino.amino_nitrogens, vec![2, 3]);
    }

    #[test]
    fn dicarbonyl_condensation_topology_keeps_one_carbonyl_explicit() {
        let center = participant_for(
            "destroy:graph:atoms=C.O.C.C.O.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-d-4,0-s-5,2-s-6,2-s-7,3-s-8",
            ReactiveSiteKind::DicarbonylElectrophile,
        )
        .dicarbonyl_electrophile_center()
        .unwrap();
        let topologies = center.condensation_topologies();
        assert_eq!(topologies.len(), 2);
        assert!(topologies.iter().any(|topology| {
            topology.retained_carbonyl_carbon == 0
                && topology.retained_carbonyl_oxygen == 1
                && topology.imine_carbon == 3
                && topology.imine_oxygen == 4
                && topology.bridge_carbon == 2
        }));
        assert!(topologies.iter().all(|topology| {
            !topology.bridge_hydrogens.is_empty() && !topology.imine_carbon_hydrogens.is_empty()
        }));
    }

    #[test]
    fn one_four_dicarbonyl_is_not_a_pyrimidine_condensation_topology() {
        let center = participant_for(
            "destroy:graph:atoms=C.O.C.C.C.O.H.H.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-s-4,4-d-5,0-s-6,2-s-7,2-s-8,3-s-9,3-s-10,4-s-11",
            ReactiveSiteKind::DicarbonylElectrophile,
        )
        .dicarbonyl_electrophile_center()
        .unwrap();
        assert!(center.condensation_topologies().is_empty());
    }
}
