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

#[derive(Clone)]
pub(crate) struct CarboxylicAcidSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) hydroxyl_oxygen: usize,
    pub(crate) hydroxyl_hydrogen: usize,
}

#[derive(Clone)]
pub(crate) struct AcylChlorideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) chlorine: usize,
}

#[derive(Clone)]
pub(crate) struct AmideSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) nitrogen: usize,
    pub(crate) nitrogen_hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct AmineSite<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) hydrogens: Vec<usize>,
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
                        "F" | "Cl" | "I"
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
        let nitrogen_hydrogens = bonded_hydrogens(self.structure, nitrogen);
        Ok(AmideSite {
            participant: self,
            carbon,
            nitrogen,
            nitrogen_hydrogens,
        })
    }

    pub(crate) fn amine_site(self) -> ChemistryResult<AmineSite<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::PrimaryAmine | ReactiveSiteKind::NonTertiaryAmine
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
