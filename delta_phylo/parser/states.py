"""
CharacterState class representing a single state of a morphological character.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class CharacterState:
    """Represents a single discrete state of a morphological character.

    Attributes:
        state_id: Numeric identifier for the state (0-based index).
        label: Human-readable label for the state.
        description: Optional longer description of the state.
    """

    state_id: int
    label: str
    description: Optional[str] = field(default=None)

    def __repr__(self) -> str:
        return f"CharacterState(id={self.state_id}, label={self.label!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CharacterState):
            return NotImplemented
        return self.state_id == other.state_id and self.label == other.label

    def __hash__(self) -> int:
        return hash((self.state_id, self.label))
