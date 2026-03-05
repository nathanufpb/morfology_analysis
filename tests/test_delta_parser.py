"""Tests for the DELTA parser."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from delta_phylo.parser.delta_reader import DeltaReader
from delta_phylo.parser.characters import Character, CharacterType
from delta_phylo.parser.taxa import Taxon
from delta_phylo.parser.states import CharacterState


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def delta_dir(tmp_path: Path) -> Path:
    """Create a minimal in-memory DELTA dataset directory."""
    specs = textwrap.dedent("""\
        *NUMBER OF CHARACTERS 4
        *MAXIMUM NUMBER OF ITEMS 3
    """)
    characters = textwrap.dedent("""\
        #1. body size/
        1. small/
        2. medium/
        3. large/

        #2. wings/
        1. absent/
        2. present/

        #3. number of legs/
        1. 6/
        2. 8/
        3. more than 8/

        #4. colour/
        1. red/
        2. blue/
        3. green/
        4. other/
    """)
    items = textwrap.dedent("""\
        #1. Taxon_A/
        1,1 2,2 3,1 4,1/

        #2. Taxon_B/
        1,2 2,1 3,2 4,3/

        #3. Taxon_C/
        1,3 2,2 3,3 4,?/
    """)
    (tmp_path / "specs").write_text(specs)
    (tmp_path / "characters").write_text(characters)
    (tmp_path / "items").write_text(items)
    return tmp_path


@pytest.fixture
def reader(delta_dir: Path) -> DeltaReader:
    r = DeltaReader(delta_dir)
    r.read()
    return r


# ---------------------------------------------------------------------------
# Tests: DeltaReader
# ---------------------------------------------------------------------------


class TestDeltaReader:
    def test_read_returns_correct_counts(self, reader: DeltaReader) -> None:
        assert len(reader.characters) == 4
        assert len(reader.taxa) == 3

    def test_character_names(self, reader: DeltaReader) -> None:
        names = [c.name for c in reader.characters]
        assert "body size" in names
        assert "wings" in names

    def test_character_ids_sequential(self, reader: DeltaReader) -> None:
        ids = [c.char_id for c in reader.characters]
        assert ids == sorted(ids)

    def test_binary_character_type(self, reader: DeltaReader) -> None:
        wings = next(c for c in reader.characters if c.char_id == 2)
        assert wings.char_type == CharacterType.BINARY

    def test_multistate_character_states(self, reader: DeltaReader) -> None:
        body_size = next(c for c in reader.characters if c.char_id == 1)
        assert body_size.num_states == 3

    def test_taxon_names(self, reader: DeltaReader) -> None:
        names = [t.name for t in reader.taxa]
        assert "Taxon_A" in names
        assert "Taxon_B" in names
        assert "Taxon_C" in names

    def test_taxon_scores_single_state(self, reader: DeltaReader) -> None:
        taxon_a = next(t for t in reader.taxa if t.name == "Taxon_A")
        # DELTA state 1 → 0-based index 0
        assert taxon_a.get_score(1) == 0
        assert taxon_a.get_score(2) == 1

    def test_missing_data(self, reader: DeltaReader) -> None:
        taxon_c = next(t for t in reader.taxa if t.name == "Taxon_C")
        assert taxon_c.is_missing(4)

    def test_specs_parsed(self, reader: DeltaReader) -> None:
        assert reader.num_characters == 4
        assert reader.num_taxa == 3


class TestDeltaReaderPolymorphic:
    """Test polymorphic state parsing."""

    def test_polymorphic_score(self, tmp_path: Path) -> None:
        (tmp_path / "characters").write_text(
            "#1. colour/\n1. red/\n2. blue/\n3. green/\n"
        )
        (tmp_path / "items").write_text(
            "#1. TaxonPoly/\n1,1/2/\n"
        )
        r = DeltaReader(tmp_path)
        r.read()
        taxon = r.taxa[0]
        score = taxon.get_score(1)
        assert isinstance(score, list)
        assert len(score) == 2

    def test_polymorphic_is_polymorphic(self, tmp_path: Path) -> None:
        (tmp_path / "characters").write_text(
            "#1. colour/\n1. red/\n2. blue/\n"
        )
        (tmp_path / "items").write_text(
            "#1. TaxonPoly/\n1,1/2/\n"
        )
        r = DeltaReader(tmp_path)
        r.read()
        taxon = r.taxa[0]
        assert taxon.is_polymorphic(1)


class TestDeltaReaderOrderedCharacters:
    """Regression tests for issue #2 — *CHARACTER TYPES OM directive."""

    def test_om_directive_sets_ordered_type(self, tmp_path: Path) -> None:
        """Characters listed under OM must be parsed as CharacterType.ORDERED."""
        (tmp_path / "specs").write_text(
            "*NUMBER OF CHARACTERS 3\n"
            "*CHARACTER TYPES OM 1 3\n"
        )
        (tmp_path / "characters").write_text(
            "#1. size/\n1. small/\n2. medium/\n3. large/\n"
            "#2. wings/\n1. absent/\n2. present/\n"
            "#3. legs/\n1. few/\n2. many/\n3. lots/\n"
        )
        (tmp_path / "items").write_text(
            "#1. TaxA/\n1,1 2,2 3,1/\n"
            "#2. TaxB/\n1,3 2,1 3,3/\n"
        )
        r = DeltaReader(tmp_path)
        r.read()
        c1 = next(c for c in r.characters if c.char_id == 1)
        c2 = next(c for c in r.characters if c.char_id == 2)
        c3 = next(c for c in r.characters if c.char_id == 3)
        assert c1.char_type == CharacterType.ORDERED, "char 1 should be ORDERED via OM"
        assert c2.char_type != CharacterType.ORDERED, "char 2 not listed in OM"
        assert c3.char_type == CharacterType.ORDERED, "char 3 should be ORDERED via OM"

    def test_ru_directive_still_works(self, tmp_path: Path) -> None:
        """Legacy RU token should also be recognised."""
        (tmp_path / "specs").write_text(
            "*NUMBER OF CHARACTERS 2\n"
            "*CHARACTER TYPES RU 2\n"
        )
        (tmp_path / "characters").write_text(
            "#1. wings/\n1. absent/\n2. present/\n"
            "#2. legs/\n1. few/\n2. many/\n3. lots/\n"
        )
        (tmp_path / "items").write_text(
            "#1. TaxA/\n1,1 2,1/\n"
            "#2. TaxB/\n1,2 2,3/\n"
        )
        r = DeltaReader(tmp_path)
        r.read()
        c2 = next(c for c in r.characters if c.char_id == 2)
        assert c2.char_type == CharacterType.ORDERED

    def test_no_character_types_directive_defaults_unordered(self, tmp_path: Path) -> None:
        (tmp_path / "characters").write_text(
            "#1. size/\n1. small/\n2. large/\n3. huge/\n"
        )
        (tmp_path / "items").write_text("#1. TaxA/\n1,1/\n")
        r = DeltaReader(tmp_path)
        r.read()
        assert r.characters[0].char_type == CharacterType.UNORDERED


class TestDeltaReaderComments:
    """Test that DELTA comments are stripped correctly."""

    def test_strip_angle_bracket_comments(self, tmp_path: Path) -> None:
        (tmp_path / "characters").write_text(
            "#1. body size <this is a comment>/\n"
            "1. small <another comment>/\n"
            "2. large/\n"
        )
        (tmp_path / "items").write_text(
            "#1. TaxonX/\n1,1/\n"
        )
        r = DeltaReader(tmp_path)
        r.read()
        char = r.characters[0]
        assert "<" not in char.name
        assert ">" not in char.name


# ---------------------------------------------------------------------------
# Tests: Character and CharacterState
# ---------------------------------------------------------------------------


class TestCharacter:
    def test_add_state(self) -> None:
        char = Character(char_id=1, name="test")
        state = char.add_state("small")
        assert state.state_id == 0
        assert state.label == "small"
        assert char.num_states == 1

    def test_get_state_by_id(self) -> None:
        char = Character(char_id=1, name="test")
        char.add_state("small")
        char.add_state("large")
        s = char.get_state_by_id(1)
        assert s is not None
        assert s.label == "large"

    def test_get_state_by_id_missing(self) -> None:
        char = Character(char_id=1, name="test")
        assert char.get_state_by_id(99) is None


# ---------------------------------------------------------------------------
# Tests: Taxon
# ---------------------------------------------------------------------------


class TestTaxon:
    def test_set_and_get_score(self) -> None:
        t = Taxon(taxon_id=1, name="TestTaxon")
        t.set_score(1, 2)
        assert t.get_score(1) == 2

    def test_missing_returns_none(self) -> None:
        t = Taxon(taxon_id=1, name="TestTaxon")
        assert t.get_score(99) is None
        assert t.is_missing(99)

    def test_character_ids(self) -> None:
        t = Taxon(taxon_id=1, name="TestTaxon")
        t.set_score(1, 0)
        t.set_score(3, 2)
        assert t.character_ids() == {1, 3}
