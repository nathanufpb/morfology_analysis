# DELTA Format Reference

DELTA (DEscription Language for TAxonomy) is a standard format for
taxonomic character descriptions, developed by M.J. Dallwitz and colleagues.

## Files

A DELTA dataset consists of three primary files (typically in a directory):

| File | Purpose |
|------|---------|
| `specs` | Dataset specifications (number of characters, taxa) |
| `characters` | Character and state definitions |
| `items` | Taxon descriptions (character scores) |

---

## specs

The `specs` file contains dataset-level directives.

```
*NUMBER OF CHARACTERS 20
*MAXIMUM NUMBER OF ITEMS 10
*CHARACTER TYPES RU 3 7 12
```

- `*NUMBER OF CHARACTERS` – total character count
- `*MAXIMUM NUMBER OF ITEMS` – total taxon count
- `*CHARACTER TYPES RU` – ordered (RU = range unit) characters

---

## characters

Character definitions use numbered entries starting with `#`.

```
#1. body size/
1. small/
2. medium/
3. large/

#2. wings/
1. absent/
2. present/
```

- Each character begins with `#<number>.`
- Character name follows the period.
- States are numbered `1.`, `2.`, etc.
- Entries end with `/`.
- Comments can be included in `<angle brackets>`.

---

## items

Taxon descriptions score each character.

```
#1. Apis mellifera <honey bee>/
1,2 2,1 5,3 7,2 10,? 12,1/2/
```

- Each taxon begins with `#<number>.`
- Taxon name follows the period.
- Scores are `<char_id>,<state>` pairs.
- States are **1-based** in DELTA files; delta-phylo converts them to **0-based** internally.
- `?` or `U` = missing data.
- `-` = inapplicable.
- `1/2` = polymorphic (both states observed).
- `1-3` = range (states 1 through 3).

---

## Character Types

| Type | Description | Example |
|------|-------------|---------|
| Unordered | Standard multistate | body colour |
| Ordered (RU) | Linear ordered series | size (small → large) |
| Binary | Two states only | presence/absence |

---

## Internal Representation

delta-phylo converts DELTA data to 0-based integer states:

| DELTA state | Internal index |
|-------------|---------------|
| 1 | 0 |
| 2 | 1 |
| 3 | 2 |
| ? | None / NaN |

---

## References

- Dallwitz, M.J., Paine, T.A., and Zurcher, E.J. (1999+). User's guide to the
  DELTA System. https://www.delta-intkey.com
- Dallwitz, M.J. (1980). A general system for coding taxonomic descriptions.
  *Taxon* 29: 41–46.
