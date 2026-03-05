"""
DeltaReader: parser for DELTA (DEscription Language for TAxonomy) datasets.

Supports reading the standard DELTA files:
  - ``characters``  (character definitions and states)
  - ``items``       (taxon descriptions / scored characters)
  - ``specs``       (dataset specifications)

Reference: Dallwitz, M.J., Paine, T.A., and Zurcher, E.J. (1999+).
  User's guide to the DELTA System. https://www.delta-intkey.com
"""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delta_phylo.parser.characters import Character, CharacterType
from delta_phylo.parser.states import CharacterState
from delta_phylo.parser.taxa import Taxon, ScoreValue

logger = logging.getLogger(__name__)


class DeltaReader:
    """Parse a directory of DELTA files into structured Python objects.

    The reader expects a directory containing at least an ``items`` file and
    (optionally) a ``characters`` and/or ``specs`` file.

    Usage::

        reader = DeltaReader("path/to/delta/dataset")
        reader.read()
        characters = reader.characters   # list[Character]
        taxa       = reader.taxa         # list[Taxon]

    Attributes:
        dataset_dir: Path to the directory containing DELTA files.
        characters: Parsed list of :class:`~delta_phylo.parser.characters.Character`.
        taxa: Parsed list of :class:`~delta_phylo.parser.taxa.Taxon`.
        num_characters: Number of characters declared in *specs*.
        num_taxa: Number of taxa declared in *specs*.
    """

    _MISSING_TOKENS: frozenset = frozenset({"?", "U", "-"})
    _INAPPLICABLE_TOKEN: str = "-"

    def __init__(self, dataset_dir: str | Path) -> None:
        self.dataset_dir = Path(dataset_dir)
        self.characters: List[Character] = []
        self.taxa: List[Taxon] = []
        self.num_characters: int = 0
        self.num_taxa: int = 0
        self._char_index: Dict[int, Character] = {}
        self._ordered_char_ids: set = set()  # populated by _read_specs

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def read(self) -> None:
        """Parse all DELTA files found in *dataset_dir*.

        Raises:
            FileNotFoundError: If no recognised DELTA files are present.
        """
        logger.info("Reading DELTA dataset from %s", self.dataset_dir)
        self._read_specs()
        self._read_characters()
        self._read_items()
        logger.info(
            "Finished reading: %d characters, %d taxa",
            len(self.characters),
            len(self.taxa),
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _find_file(self, *names: str) -> Optional[Path]:
        """Return the first existing file matching any of *names* (case-insensitive)."""
        for name in names:
            for entry in self.dataset_dir.iterdir():
                if entry.is_file() and entry.name.lower() == name.lower():
                    return entry
        return None

    @staticmethod
    def _load_text(path: Path) -> str:
        """Read a file, trying UTF-8 then latin-1 encoding."""
        for enc in ("utf-8", "latin-1"):
            try:
                return path.read_text(encoding=enc)
            except UnicodeDecodeError:
                continue
        return path.read_text(errors="replace")

    @staticmethod
    def _strip_comments(text: str) -> str:
        """Remove DELTA comments enclosed in angle brackets ``<...>``."""
        # Remove nested angle-bracket comments greedily.
        # DELTA allows nesting; we iterate until stable.
        prev = None
        while prev != text:
            prev = text
            text = re.sub(r"<[^<>]*>", "", text)
        return text

    @staticmethod
    def _normalise(text: str) -> str:
        """Normalise whitespace and remove trailing/leading blanks."""
        return " ".join(text.split())

    @staticmethod
    def _parse_ordered_ids(text: str) -> set:
        """Extract character IDs marked as ordered from a *CHARACTER TYPES block.

        Supports both the canonical ``OM`` (Ordered Multistate) token and the
        legacy ``RU`` token used by some older DELTA tools.

        Args:
            text: Raw (comment-stripped) text to search within.

        Returns:
            Set of 1-based character IDs that are ordered.
        """
        ordered: set = set()
        types_match = re.search(
            r"\*CHARACTER\s+TYPES\s*(.*?)(?=\*|\Z)", text, re.IGNORECASE | re.DOTALL
        )
        if not types_match:
            return ordered
        for m in re.finditer(
            r"(?:OM|RU)\s+([\d,/\s]+)", types_match.group(1), re.IGNORECASE
        ):
            for part in re.split(r"[,\s]+", m.group(1)):
                part = part.strip()
                if part.isdigit():
                    ordered.add(int(part))
        return ordered

    # ------------------------------------------------------------------
    # specs
    # ------------------------------------------------------------------

    def _read_specs(self) -> None:
        """Parse the ``specs`` file for dataset-level metadata."""
        path = self._find_file("specs")
        if path is None:
            logger.debug("No 'specs' file found; skipping.")
            return
        logger.debug("Parsing specs: %s", path)
        text = self._strip_comments(self._load_text(path))
        # *NC directive: number of characters
        m = re.search(r"\*NUMBER\s+OF\s+CHARACTERS\s+(\d+)", text, re.IGNORECASE)
        if m:
            self.num_characters = int(m.group(1))
        # *NT directive: number of taxa (items)
        m = re.search(r"\*MAXIMUM\s+NUMBER\s+OF\s+ITEMS\s+(\d+)", text, re.IGNORECASE)
        if m:
            self.num_taxa = int(m.group(1))
        # *CHARACTER TYPES OM/RU <ids>  — ordered multistate characters
        self._ordered_char_ids = self._parse_ordered_ids(text)
        logger.debug(
            "specs: num_characters=%d, num_taxa=%d, ordered=%s",
            self.num_characters,
            self.num_taxa,
            self._ordered_char_ids,
        )

    # ------------------------------------------------------------------
    # characters
    # ------------------------------------------------------------------

    def _read_characters(self) -> None:
        """Parse the ``characters`` file into :class:`Character` objects."""
        path = self._find_file("characters", "chars")
        if path is None:
            logger.debug("No 'characters' file found; creating placeholders.")
            # Create placeholder characters based on specs count
            for i in range(1, self.num_characters + 1):
                char = Character(char_id=i, name=f"Character {i}")
                self.characters.append(char)
                self._char_index[i] = char
            return
        logger.debug("Parsing characters: %s", path)
        raw = self._load_text(path)
        text = self._strip_comments(raw)
        self._parse_character_list(text)

    def _parse_character_list(self, text: str) -> None:
        """Extract character definitions from the character list text."""
        # Split on character markers like  #1.  or  #1/
        # DELTA uses  #<number>.  to start character definitions.
        char_blocks = re.split(r"(?=#\d+\.)", text)

        # Ordered character IDs come from *CHARACTER TYPES in specs (parsed
        # during _read_specs) or, as a fallback, from the characters file itself.
        ordered_chars: set = self._ordered_char_ids.copy()
        if not ordered_chars:
            ordered_chars = self._parse_ordered_ids(text)

        for block in char_blocks:
            block = block.strip()
            if not block:
                continue
            # Match  #<id>.  at the start
            m = re.match(r"#(\d+)\.\s*(.*)", block, re.DOTALL)
            if not m:
                continue
            char_id = int(m.group(1))
            rest = m.group(2).strip()

            # Split name from states (states start with a numbered list or /)
            # States are listed as:  1. <state>  2. <state> ...
            parts = re.split(r"\b(\d+)\.\s+", rest)
            # parts[0] = character name, parts[1::2] = state numbers, parts[2::2] = labels
            char_name = self._normalise(parts[0].strip().strip("/").strip())
            if not char_name:
                char_name = f"Character {char_id}"

            char_type = (
                CharacterType.ORDERED
                if char_id in ordered_chars
                else CharacterType.UNORDERED
            )
            char = Character(char_id=char_id, name=char_name, char_type=char_type)

            # Parse states
            state_labels = parts[2::2]  # every second element starting from index 2
            for i, label in enumerate(state_labels):
                label_clean = self._normalise(label.strip("/").strip())
                if label_clean:
                    char.add_state(label=label_clean)

            # If only 2 states, mark as BINARY
            if char.num_states == 2 and char_type == CharacterType.UNORDERED:
                char.char_type = CharacterType.BINARY

            self.characters.append(char)
            self._char_index[char_id] = char
            logger.debug("Parsed %r", char)

        # Fill in any missing characters up to num_characters
        existing_ids = {c.char_id for c in self.characters}
        for i in range(1, self.num_characters + 1):
            if i not in existing_ids:
                char = Character(char_id=i, name=f"Character {i}")
                self.characters.append(char)
                self._char_index[i] = char

        # Sort by char_id
        self.characters.sort(key=lambda c: c.char_id)

    # ------------------------------------------------------------------
    # items (taxa)
    # ------------------------------------------------------------------

    def _read_items(self) -> None:
        """Parse the ``items`` file into :class:`Taxon` objects."""
        path = self._find_file("items")
        if path is None:
            logger.warning("No 'items' file found in %s", self.dataset_dir)
            return
        logger.debug("Parsing items: %s", path)
        raw = self._load_text(path)
        text = self._strip_comments(raw)
        self._parse_item_list(text)

    def _parse_item_list(self, text: str) -> None:
        """Extract taxon descriptions from the items text."""
        # Items start with  #<number>.  (same as characters format)
        item_blocks = re.split(r"(?=#\d+\.)", text)
        for block in item_blocks:
            block = block.strip()
            if not block:
                continue
            m = re.match(r"#(\d+)\.\s*(.*)", block, re.DOTALL)
            if not m:
                continue
            taxon_id = int(m.group(1))
            rest = m.group(2).strip()

            # Taxon name is first text before the first character score.
            # Character scores are given as  <char_id>,<score>/
            # or  <char_id>,<score> <char_id2>,<score2> ...
            # The name ends at the first  digit,  pattern (score entry).
            name_match = re.match(r"^(.*?)(?=\d+,)", rest, re.DOTALL)
            if name_match:
                taxon_name = self._normalise(
                    name_match.group(1).strip().strip("/").strip()
                )
                scores_text = rest[name_match.end():]
            else:
                taxon_name = self._normalise(rest)
                scores_text = ""

            if not taxon_name:
                taxon_name = f"Taxon_{taxon_id}"

            taxon = Taxon(taxon_id=taxon_id, name=taxon_name)

            # Parse character scores
            self._parse_scores(taxon, scores_text)
            self.taxa.append(taxon)
            logger.debug("Parsed %r with %d scores", taxon, len(taxon.scores))

        # Sort taxa by ID
        self.taxa.sort(key=lambda t: t.taxon_id)

    def _parse_scores(self, taxon: Taxon, scores_text: str) -> None:
        """Parse individual character-score entries into a Taxon's scores dict.

        DELTA score syntax examples:
          ``1,0``          single state 0 for character 1
          ``3,1/2``        polymorphic: states 1 and 2 for character 3
          ``5,?``          missing data for character 5
          ``7,-``          inapplicable for character 7
          ``2,0-2``        states 0 through 2 (range; treated as polymorphic)

        Args:
            taxon: Taxon object to populate.
            scores_text: Raw string of scored characters.
        """
        # Match patterns like  <char_id>,<value>
        pattern = re.compile(
            r"(\d+),\s*"           # character id
            r"([\d\.\-/\?U]+)"     # score value(s)
        )
        for m in pattern.finditer(scores_text):
            char_id = int(m.group(1))
            raw_val = m.group(2).strip()
            value = self._decode_score(raw_val)
            taxon.set_score(char_id, value)

    @staticmethod
    def _decode_score(raw: str) -> ScoreValue:
        """Convert a raw DELTA score string to a Python value.

        Args:
            raw: Raw score string from DELTA file.

        Returns:
            ``None`` for missing/inapplicable data,
            ``int`` for a single state,
            ``list[int]`` for polymorphic states.
        """
        raw = raw.strip()
        # Missing / inapplicable
        if raw in ("?", "U", "-"):
            return None
        # Polymorphic: slash-separated  e.g.  0/1  or  1/2/3
        if "/" in raw:
            parts = [p.strip() for p in raw.split("/") if p.strip()]
            states: List[int] = []
            for part in parts:
                if part.isdigit():
                    # Convert from 1-based DELTA state numbers to 0-based index
                    states.append(int(part) - 1)
            return states if len(states) > 1 else (states[0] if states else None)
        # Range: e.g.  1-3  → states 1,2,3
        range_m = re.match(r"^(\d+)-(\d+)$", raw)
        if range_m:
            lo, hi = int(range_m.group(1)), int(range_m.group(2))
            # Convert from 1-based to 0-based
            return list(range(lo - 1, hi))
        # Single integer state (1-based in DELTA → 0-based internally)
        if raw.isdigit():
            return int(raw) - 1
        return None
