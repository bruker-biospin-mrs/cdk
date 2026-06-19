# CDK 2.12 Upgrade Notes — Branch `cdk-2.12-BBIO`

This document records all merge conflict resolutions and post-merge fixes made when merging
`cdk-2.12` (upstream CDK 2.12) into the BBIO-patched branch (`origin/ModularV2000`).

The resulting branch is `cdk-2.12-BBIO` (version `2.12-BBIO-SNAPSHOT`).

---

## Overview

The BBIO patch set refactors `MDLV2000Reader` to delegate to pluggable handler classes:

- `V2000MoleculeBlockHandler` / `V2000SlowMoleculeBlockHandler`
- `V2000PropertiesBlockHandler` / `V2000SlowPropertiesBlockHandler`
- `V2000NonStructuralDataBlockHandler`

These are exposed via factory methods on `MDLV2000Reader`:
- `newMoleculeBlockHandler()`
- `newPropertiesBlockHandler()`
- `newNonStructuralDataBlockHandler()`

`IteratingSDFReader` was extended with a `Builder` supporting custom reader suppliers,
end-of-molecule detection, and reading mode.

---

## Merge: `cdk-2.12` → `cdk-2.12-BBIO`

Base commit: merge commit `28576ad54f` (merge of `cdk-2.12` into `ModularV2000`).

### Files with conflicts (all in `storage/ctab`)

#### `MDLV2000Reader.java`
**Resolution:** Started from CDK 2.12, then applied BBIO additions on top.

BBIO additions kept:
- `getReaderMode()` — exposes the reader mode to handler classes
- `newMoleculeBlockHandler()` — factory method for overriding molecule block handling
- `newPropertiesBlockHandler()` — factory method for overriding properties block handling
- `newNonStructuralDataBlockHandler()` — factory method for overriding non-structural data handling
- The `ATOM_ALIAS` case in `readPropertiesFast()` was wired to call `propHandler.handleAtomAlias()`

CDK 2.12 changes kept:
- All stereo perception improvements
- `IChemObject.AROMATIC` usage (replacing `CDKConstants.ISAROMATIC` in molecule block)
- Improved error handling and lenient parsing

#### `MDLV3000Reader.java`
**Resolution:** Took CDK 2.12 entirely (`--theirs`). No BBIO-specific changes were needed
because:
- The SGroup crash (MFSW-1822) is fixed in CDK 2.12 (`readSGroup` is now private and correct)
- CDK 2.12 already uses `IChemObject.AROMATIC` for aromaticity

#### `iterator/IteratingSDFReader.java`
**Resolution:** Started from CDK 2.12, then applied BBIO additions on top.

CDK 2.12's `hasNext()` was rewritten compared to the BBIO version. The BBIO approach — reading
until `endOfMoleculeFunction` returns true, then calling `readDataBlockInto()` for SDF properties
— was retained because CDK 2.12's `hasNext()` did not read SDF data properties at all.

BBIO additions kept:
- `Builder` inner class with `setV2000ReaderSupplier()`, `setV3000ReaderSupplier()`,
  `setEndOfMoleculeFunction()`, `setSkip()`, `setReadingMode()`
- `readDataBlockInto()` and `extractFieldData()` for SDF property reading
- `getReaderMode()` accessor

Uses `java.util.function.Supplier` and `java.util.function.Function` (not Guava).

#### Test files (6 files)
All test files were taken from CDK 2.12 (`--theirs`):
- `MDLV2000ReaderTest.java`
- `MDLV3000ReaderTest.java`
- `iterator/IteratingSDFReaderTest.java`
- `MDLV2000AtomBlockTest.java`
- `MDLV2000BondBlockTest.java`
- `MDLV2000PropertiesBlockTest.java`

---

## Post-Merge Bug Fixes

Three test failures were found after the merge and fixed:

### Fix 1 — `V2000PropertiesBlockHandler.label()`: preserve atomic number
**File:** `storage/ctab/src/main/java/org/openscience/cdk/io/V2000PropertiesBlockHandler.java`

**Test:** `MDLV2000ReaderTest.keepAtomicNumberOfAlias`

**Problem:** When replacing a real atom with a pseudo atom for an atom alias,
`pseudoAtom.setAtomicNumber(0)` discarded the original atom's atomic number.
CDK 2.12's `MDLV2000Reader.label()` preserves it.

**Fix:** Changed `pseudoAtom.setAtomicNumber(0)` to `pseudoAtom.setAtomicNumber(atom.getAtomicNumber())`.

---

### Fix 2 — `V2000BlockHandler.readMolfileInt()`: bounds check for short lines
**File:** `storage/ctab/src/main/java/org/openscience/cdk/io/V2000BlockHandler.java`

**Test:** `MDLV2000ReaderTest.test` (line 1844)

**Problem:** The atom alias format `"A  aaa"` normally has 6+ characters. However, a test mol
file contained `"A   1"` (5 chars). `readMolfileInt(line, 3)` accesses `line[3]`, `line[4]`,
`line[5]` — but index 5 is out of bounds for a 5-char string.

CDK 2.12's `MDLV2000Reader.readMolfileInt()` guards each access:
```java
if (index + 1 == line.length()) return sign * result;
if (index + 2 == line.length()) return sign * result;
```

**Fix:** Added the same two bounds checks to `V2000BlockHandler.readMolfileInt()`.

---

### Fix 3 — `IteratingSDFReader.extractFieldData()`: trim all trailing blank lines
**File:** `storage/ctab/src/main/java/org/openscience/cdk/io/iterator/IteratingSDFReader.java`

**Test:** `IteratingSDFReaderTest.testExtraSpaces`

**Problem:** The ChEMBL API test SDF file (`chemblApiExamples.sdf`) has two consecutive blank
lines after each property value, rather than the standard single blank line. This caused
`extractFieldData()` to build a value like `"MILCICLIB\n\n"`. The original single-trim
```java
if (len > 1 && data.charAt(len - 1) == '\n') data.setLength(len - 1);
```
only removed one trailing newline, leaving `"MILCICLIB\n"` instead of `"MILCICLIB"`.

**Fix:** Changed to a loop:
```java
while (data.length() > 0 && data.charAt(data.length() - 1) == '\n')
    data.setLength(data.length() - 1);
```

---

## `V2000PropertiesBlockHandler` — CDK 2.12 New Property Cases

The following property line types were added to `V2000PropertiesBlockHandler.readProperties()`
to match CDK 2.12 behaviour (they existed in `MDLV2000Reader` but were not in the BBIO handler):

| Case          | Handler method added           | Notes                                       |
|---------------|-------------------------------|---------------------------------------------|
| `LEGACY_ATOM_LIST` | `handleLegacyAtomList()`  | `A  aaa` atom list (old V2000 format)       |
| `M_ALS`       | `handleAtomListSgroup()`       | M ALS — extended atom list                  |
| `M_SDT`       | `handleSgroupDataType()`       | Sgroup data type                            |
| `M_SDD`       | `handleSgroupDataDisplay()`    | Sgroup data display                         |
| `M_SCD`/`M_SED` | `handleSgroupDataComplete()` | Sgroup data (multi-line, terminated by SED) |

---

## `V2000BlockHandler.java` — Guava Removal

`ImmutableSet` (Guava) was replaced with `Collections.unmodifiableSet(new HashSet<>(Arrays.asList(...)))` for Java 8 compatibility without Guava dependency.

---

## tsjava Changes Required

The following changes were made to `tsjava/common_lib` to compile against `2.12-BBIO-SNAPSHOT`:

### `common_lib/pom.xml` (via parent `tsjava/pom.xml`)
Changed `cdk.version` from `2.0.2-BBIO` to `2.12-BBIO-SNAPSHOT`.

### `SDFReader.java`
Removed `ConfiguredMDLV3000Reader.readSGroup(IAtomContainer)` override.

**Reason:** In CDK 2.12, `MDLV3000Reader.readSGroup` is now `private void readSGroup(ReadState)`
(private, different parameter type). The old `public void readSGroup(IAtomContainer)` no longer
exists, so `@Override` would fail to compile. The underlying SGroup crash (MFSW-1822) is fixed
in CDK 2.12, so the workaround override is no longer needed.

**Note:** `CDKConstants.ISAROMATIC` is deprecated in CDK 2.12 in favour of `IChemObject.AROMATIC`
(same int value `0x0020`). The deprecated constant is still present and compiles. It is used in
`SDFReader.kekulizeIfAromatic()` and in the BBIO handler classes — these are intentional leftovers
that can be cleaned up in a future pass.
