use crate::chemistry::error::ChemistryResult;
use crate::chemistry::reaction::Reaction;

use super::resolver::DerivedSubstanceResolver;

pub(crate) trait ReactionSolver {
    type Input;

    fn solve(
        input: Self::Input,
        resolver: &mut DerivedSubstanceResolver,
    ) -> ChemistryResult<Vec<Reaction>>;
}

pub(crate) trait PairReactionSolver {
    type Left;
    type Right;

    fn solve_pair(
        left: Self::Left,
        right: Self::Right,
        resolver: &mut DerivedSubstanceResolver,
    ) -> ChemistryResult<Vec<Reaction>>;
}
