# FreeBayes Test Suite

## Overview

This directory contains tests for FreeBayes, including:
- **Bash TAP tests** (`t/*.t`) - Integration tests for command-line behavior
- **C++ unit tests** (`test_*.cpp`) - Unit tests for core components

## Test Categories

### 1. Baseline Tests (`test_allele_parser_baseline.cpp`)

**Purpose:** Verify current AlleleParser behavior before refactoring.

These tests establish a baseline for AlleleParser functionality:
- BAM file loading and parsing
- Reference sequence access
- Alignment iteration
- Allele detection
- Region handling
- Coverage tracking
- Parameter handling

**Why:** Before refactoring AlleleParser to use the new IAlignmentReader interface,
we need tests that verify the current behavior. These tests will ensure the refactoring
doesn't break existing functionality.

**Test Data:**
- `data/test.bam` - Minimal test with 1 read (ref: "ATCGGCTAAAA", 11bp)
- `data/test.ref` - Reference sequence
- `tiny/NA12878.chr22.tiny.bam` - Real reads from chr22
- `tiny/q.fa` - Reference for tiny BAM

### 2. Interface Tests (`test_alignment_interfaces.cpp`)

**Purpose:** Verify the new alignment reader abstraction layer.

Tests for the decoupled alignment reading interfaces:
- Format detection (by extension and magic bytes)
- Reader creation via factory
- Interface implementation (SeqLibAlignmentReader)
- Polymorphic usage

**Why:** Ensures the new interface layer works correctly before integrating
with AlleleParser.

### 3. Integration Tests (Bash TAP tests in `t/`)

Command-line integration tests:
- `00_region_and_target_handling.t` - Region/target specification
- `01_call_variants.t` - Basic variant calling
- `02_multi_bam.t` - Multiple BAM file handling
- `03_reference_bases.t` - Reference base handling

## Running Tests

### All Tests (via Meson)

```bash
# From freebayes root directory
cd build/
ninja
meson test
```

### Specific Test Suites

```bash
# Run only C++ unit tests
meson test test_allele_parser_baseline
meson test test_alignment_interfaces

# Run only bash TAP tests
meson test T00
meson test T01
```

### Individual Test Executables

```bash
# Build and run baseline tests
cd build/
ninja test_allele_parser_baseline
./test_allele_parser_baseline

# Build and run interface tests
ninja test_alignment_interfaces
./test_alignment_interfaces
```

### Verbose Output

```bash
# See detailed test output
meson test --verbose test_allele_parser_baseline
```

## Test Data Structure

```
test/
├── data/                          # Minimal test data
│   ├── test.bam                  # 1 read with SNPs
│   ├── test.bam.bai              # Index
│   ├── test.ref                  # 11bp reference
│   ├── test.bed                  # Target regions
│   └── test.vcf                  # Expected variants
│
├── tiny/                          # Small real-world data
│   ├── NA12878.chr22.tiny.bam    # Chr22 reads
│   ├── NA12878.chr22.tiny.cram   # CRAM version
│   ├── q.fa                       # Reference
│   └── *.giab.vcf                # GIAB truth set
│
└── t/                             # Bash TAP tests
    ├── 00_region_and_target_handling.t
    ├── 01_call_variants.t
    └── ...
```

## Adding New Tests

### C++ Unit Test

1. Create `test_mytest.cpp` in `test/` directory
2. Use the test framework macros:
   ```cpp
   TEST(MyTestName) {
       // Test code
       ASSERT_EQ(actual, expected);
   }

   int main() {
       RUN_TEST(MyTestName);
       return 0;
   }
   ```

3. Add to `meson.build`:
   ```meson
   test_mytest = executable('test_mytest', ...)
   test('test_mytest', test_mytest, workdir: testdir)
   ```

### Bash TAP Test

1. Create `t/NN_description.t`
2. Use bash-tap framework:
   ```bash
   #!/usr/bin/env bash
   . bash-tap/bash-tap-bootstrap

   plan tests 3

   is $(freebayes ...) "expected" "test description"
   ```

3. Add to `meson.build`:
   ```meson
   test('TNN', prove, args: ['-e','bash','-v','t/NN_description.t'])
   ```

## Test-Driven Refactoring Plan

### Phase 1: Establish Baseline ✅
- [x] Create `test_allele_parser_baseline.cpp`
- [x] Verify all tests pass with current implementation
- [x] Tests cover: BAM reading, allele detection, region handling

### Phase 2: Interface Implementation ✅
- [x] Create interface layer (IAlignment, IAlignmentReader)
- [x] Create `test_alignment_interfaces.cpp`
- [x] Verify interfaces work in isolation

### Phase 3: Integration (TODO)
- [ ] Add interface member to AlleleParser
- [ ] Replace BAM-specific calls with interface calls
- [ ] Verify `test_allele_parser_baseline.cpp` still passes
- [ ] No changes to test code needed (tests interface, not implementation)

### Phase 4: Extended Testing (TODO)
- [ ] Add tests for new format support (PAF, etc.)
- [ ] Add tests for mock readers
- [ ] Performance regression tests

## Expected Test Results

### Current Status

With the **existing AlleleParser** implementation:
- ✅ `test_allele_parser_baseline` should pass
- ✅ `test_alignment_interfaces` should pass (interfaces only)
- ✅ All bash TAP tests (`T00-T03`) should pass

### After Refactoring

With **refactored AlleleParser** (using IAlignmentReader):
- ✅ `test_allele_parser_baseline` should still pass (same behavior)
- ✅ `test_alignment_interfaces` should still pass
- ✅ All bash TAP tests should still pass (same command-line behavior)
- ✅ New tests for format flexibility (e.g., mock readers)

## Debugging Failed Tests

### Test Fails to Compile

```bash
# Check compilation errors
cd build/
ninja test_allele_parser_baseline 2>&1 | less
```

Common issues:
- Missing includes
- Linker errors (check dependencies in meson.build)
- API changes in AlleleParser

### Test Fails at Runtime

```bash
# Run with verbose output
./test_allele_parser_baseline

# Or with meson
meson test --verbose test_allele_parser_baseline
```

Common issues:
- Test data not found (check working directory)
- BAM index missing (regenerate with `samtools index`)
- Assertions failing (behavior changed)

### Test Data Issues

```bash
# Verify test BAM is valid
samtools view test/data/test.bam | head

# Check reference
cat test/data/test.ref

# Regenerate index if needed
samtools index test/data/test.bam
```

## Test Coverage

Current test coverage for AlleleParser:

| Component | Coverage | Notes |
|-----------|----------|-------|
| BAM Reading | ✅ High | Multiple tests |
| Reference Access | ✅ High | Tests with known reference |
| Allele Detection | ✅ Medium | Tests with real data |
| Region Handling | ✅ High | Multiple region tests |
| Parameter Parsing | ✅ Medium | Tests key parameters |
| Coverage Tracking | ✅ Low | Basic test only |
| Genotyping | ⚠️ None | Future work |
| VCF Output | ⚠️ None | Future work |

## Continuous Integration

Tests should be run:
- Before committing changes
- In CI/CD pipeline
- Before releases

Recommended CI setup:
```yaml
- name: Build
  run: meson setup build/ && cd build && ninja

- name: Run Tests
  run: cd build && meson test --verbose
```

## Performance Testing

For performance regression testing:
```bash
# Time variant calling
time freebayes -f test/tiny/q.fa test/tiny/NA12878.chr22.tiny.bam > /dev/null

# Compare before/after refactoring
# Should be within 5% of baseline
```

## Further Reading

- [Bash TAP Framework](https://github.com/illusori/bash-tap)
- [Meson Test Guide](https://mesonbuild.com/Unit-tests.html)
- Architecture: `docs/ALIGNMENT_ABSTRACTION_DESIGN.md`
- Implementation: `docs/INTERFACE_IMPLEMENTATION_SUMMARY.md`
