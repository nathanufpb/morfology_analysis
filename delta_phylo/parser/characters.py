"""
Character and CharacterType classes for morphological characters.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Optional

from delta_phylo.parser.states import CharacterState


class CharacterType(Enum):
    """Enumeration of morphological character types."""

    UNORDERED = auto()   # Standard multistate, unordered (default)
    ORDERED = auto()     # Ordered multistate (linear series)
    BINARY = auto()      # Two-state character (0/1)
    NUMERIC = auto()     # Continuous/numeric character


@dataclass
class Character:
    """Represents a morphological character from a DELTA dataset.

    Attributes:
        char_id: 1-based identifier matching the DELTA file numbering.
        name: Full text description of the character.
        char_type: Type of the character (ordered, unordered, binary).
        states: List of possible states for this character.
        weight: Analytical weight (default 1.0).
        is_applicable: Whether the character is generally applicable.
    """

    char_id: int
    name: str
    char_type: CharacterType = CharacterType.UNORDERED
    states: List[CharacterState] = field(default_factory=list)
    weight: float = 1.0
    is_applicable: bool = True

    def __repr__(self) -> str:
        return (
            f"Character(id={self.char_id}, name={self.name!r}, "
            f"type={self.char_type.name}, states={len(self.states)})"
        )

    @property
    def num_states(self) -> int:
        """Return the number of defined states."""
        return len(self.states)

    def get_state_by_id(self, state_id: int) -> Optional[CharacterState]:
        """Retrieve a state by its numeric ID.

        Args:
            state_id: 0-based state identifier.

        Returns:
            CharacterState if found, else None.
        """
        for state in self.states:
            if state.state_id == state_id:
                return state
        return None

    def add_state(self, label: str, description: Optional[str] = None) -> CharacterState:
        """Add a new state to this character.

        Args:
            label: Short label for the state.
            description: Optional longer description.

        Returns:
            The newly created CharacterState.
        """
        state_id = len(self.states)
        state = CharacterState(state_id=state_id, label=label, description=description)
        self.states.append(state)
        return state
