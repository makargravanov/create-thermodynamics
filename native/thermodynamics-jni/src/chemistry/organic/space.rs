use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularStructure;
use crate::chemistry::reactive_site::{try_find_reactive_sites, ReactiveSite, ReactiveSiteKind};
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::{Substance, SubstanceId};

pub(crate) struct GenerationScope {
    substances: BTreeSet<SubstanceId>,
}

impl GenerationScope {
    pub(crate) fn all(registry: &ChemistryRegistry) -> Self {
        Self {
            substances: registry
                .substances()
                .map(|substance| substance.id.clone())
                .collect(),
        }
    }

    #[cfg(test)]
    pub(crate) fn from_substances(substances: &BTreeSet<SubstanceId>) -> Self {
        Self {
            substances: substances.clone(),
        }
    }

    fn contains(&self, substance_id: &SubstanceId) -> bool {
        self.substances.contains(substance_id)
    }
}

#[derive(Clone)]
pub(crate) struct SiteParticipant<'a> {
    pub(crate) substance: &'a Substance,
    pub(crate) structure: &'a MolecularStructure,
    pub(crate) site: ReactiveSite,
}

impl<'a> SiteParticipant<'a> {
    pub(crate) fn is_seed(&self, seeds: Option<&BTreeSet<SubstanceId>>) -> bool {
        seeds.is_none_or(|seeds| seeds.contains(&self.substance.id))
    }
}

pub(crate) struct OrganicGenerationSpace<'a> {
    pub(crate) all_substances: Vec<&'a Substance>,
    scope_substances: BTreeSet<SubstanceId>,
    participants_by_site: BTreeMap<ReactiveSiteKind, Vec<SiteParticipant<'a>>>,
}

impl<'a> OrganicGenerationSpace<'a> {
    pub(crate) fn new(
        substances: impl IntoIterator<Item = &'a Substance>,
        scope: &GenerationScope,
    ) -> ChemistryResult<Self> {
        let mut all_substances = Vec::new();
        let mut participants_by_site: BTreeMap<ReactiveSiteKind, Vec<SiteParticipant<'a>>> =
            BTreeMap::new();

        for substance in substances {
            all_substances.push(substance);
            if !scope.contains(&substance.id) {
                continue;
            }
            let Some(structure) = substance.molecular_structure.as_ref() else {
                continue;
            };
            for site in try_find_reactive_sites(structure)? {
                participants_by_site
                    .entry(site.kind.clone())
                    .or_default()
                    .push(SiteParticipant {
                        substance,
                        structure,
                        site,
                    });
            }
        }

        Ok(Self {
            all_substances,
            scope_substances: scope.substances.clone(),
            participants_by_site,
        })
    }

    pub(crate) fn from_substances_for_scope(
        substances: impl IntoIterator<Item = &'a Substance>,
        scope: &BTreeSet<SubstanceId>,
    ) -> ChemistryResult<Self> {
        Self::new(
            substances,
            &GenerationScope {
                substances: scope.clone(),
            },
        )
    }

    pub(crate) fn site_participants(&self) -> impl Iterator<Item = SiteParticipant<'a>> + '_ {
        self.participants_by_site
            .values()
            .flat_map(|participants| participants.iter().cloned())
    }

    pub(crate) fn sites_of(
        &self,
        kind: &ReactiveSiteKind,
    ) -> impl Iterator<Item = SiteParticipant<'a>> + '_ {
        self.participants_by_site
            .get(kind)
            .into_iter()
            .flat_map(|participants| participants.iter().cloned())
    }

    pub(crate) fn contains_substance(&self, substance_id: &str) -> bool {
        self.scope_substances
            .contains(&SubstanceId::from(substance_id))
    }
}
