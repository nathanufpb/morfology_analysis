"""
Taxon class representing a single operational taxonomic unit (OTU).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Union

# Type alias for a single character score entry.
# Values may be: int, list of ints (polymorphic), or None (missing).
ScoreValue = Optional[Union[int, List[int]]]


@dataclass
class Taxon:
    """Represents an operational taxonomic unit (OTU) in a DELTA dataset.

    Attributes:
        taxon_id: 1-based identifier matching DELTA item numbering.
        name: Scientific or common name of the taxon.
        scores: Mapping from character ID (1-based) to observed state(s).
            - ``None`` means missing data (``?``).
            - An ``int`` means a single observed state.
            - A ``list`` of ints means polymorphic states (e.g. ``{0,1}``).
    """

    taxon_id: int
    name: str
    scores: Dict[int, ScoreValue] = field(default_factory=dict)

    def __repr__(self) -> str:
        return f"Taxon(id={self.taxon_id}, name={self.name!r})"

    def get_score(self, char_id: int) -> ScoreValue:
        """Return the observed score for a given character.

        Args:
            char_id: 1-based character identifier.

        Returns:
            State value, list of states, or None for missing data.
        """
        return self.scores.get(char_id)

    def set_score(self, char_id: int, value: ScoreValue) -> None:
        """Set the observed score for a given character.

        Args:
            char_id: 1-based character identifier.
            value: State value, list of states, or None for missing.
        """
        self.scores[char_id] = value

    def is_missing(self, char_id: int) -> bool:
        """Check whether data for a character is missing.

        Args:
            char_id: 1-based character identifier.

        Returns:
            True if the score is missing (None) or not recorded.
        """
        return self.scores.get(char_id) is None

    def is_polymorphic(self, char_id: int) -> bool:
        """Check whether a taxon is polymorphic for a character.

        Args:
            char_id: 1-based character identifier.

        Returns:
            True if multiple states were observed.
        """
        val = self.scores.get(char_id)
        return isinstance(val, list) and len(val) > 1

    def character_ids(self) -> Set[int]:
        """Return the set of character IDs that have been scored.

        Returns:
            Set of character IDs.
        """
        return set(self.scores.keys())
