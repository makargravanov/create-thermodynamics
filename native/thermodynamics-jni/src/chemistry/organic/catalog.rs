use crate::chemistry::reaction::Reaction;
use crate::chemistry::substance::Substance;

#[derive(Debug, Default)]
pub(crate) struct GeneratedOrganicCatalog {
    pub(crate) substances: Vec<Substance>,
    pub(crate) reactions: Vec<Reaction>,
}
