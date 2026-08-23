# 5-point-titratiion
Convert pascal code into python

This repository contains the listing of the source code for the source code file. 

The program was coded using Turbo Pascal Ver 4.0. It allows the calculation of H2C03*alkalinity, SCFA (as AT) and systematic pH measurement error from the data collected in the 5 pH point titration procedure. (For 5 pH point titration procedure, see Appendix V).

# Repository Analysis
## 5-point-titration
### Application Type

Scientific computing library — a pure computational package with no web UI, database, or external service dependencies. The original Pascal program was a DOS text-mode interactive application (using Turbo Pascal's crt unit for screen drawing); the Python modernization strips the UI layer and exposes the calculation algorithms as a Python package with an optional CLI.
Origin and Domain Context

The source code comes from Appendix W of a 1991 doctoral thesis by Moosbrugger at the University of Cape Town (UCT). The program implements a 5-pH-point titration method for determining:

    H2CO3 alkalinity* (mg/L as CaCO3) — a measure of carbonate buffering capacity
    Short-chain fatty acid (SCFA) concentration (mg/L as acetic acid) — important for anaerobic digester monitoring
    Systematic pH measurement error — estimated via iterative convergence of two independent Ct calculations

The titration procedure measures pH at five points: initial pH (pH0) plus four target regions (~6.7, ~5.9, ~5.2, ~4.3), with known volumes of strong acid titrant added between each.
Repository Structure

5-point-titration/
├── .github/workflows/python-tests.yml   # CI: pytest on Python 3.10/3.11/3.12
├── .gitignore
├── README.md
├── pyproject.toml                        # Build config (setuptools), deps, pytest config
├── requirements.txt                      # Flat dependency list
├── config/.gitkeep                       # Placeholder for future config files
├── data/
│   ├── expected-output/.gitkeep          # Placeholder for reference outputs
│   └── sample-input/.gitkeep            # Placeholder for sample inputs
├── docs/thesis-reference/
│   ├── 5_point_titration_code.pas        # OCR'd Pascal source (primary reference)
│   ├── 1991_moosbrugger_re_5_point_titration_code.txt  # Text version of the thesis appendix
│   ├── 5point_titration.py               # Draft Python translation #1 (incomplete)
│   ├── 5point_titration_ver1.py          # Draft Python translation #2 (incomplete)
│   ├── 1991_moosbrugger_re_5_point_titration_code.docx
│   ├── 1991_moosbrugger_re_%_point_titration_code.pdf
│   ├── 5 Point Titration.docx
│   └── 5-Point.Titration.pdf
├── legacy-code/pascal/.gitkeep           # Placeholder (actual Pascal is in docs/)
├── scripts/.gitkeep                      # Placeholder for utility scripts
├── src/titration/
│   ├── __init__.py                       # Package init, exports, version 0.1.0
│   └── core.py                           # Stub implementations (placeholders only)
└── tests/
    ├── __init__.py
    ├── conftest.py                       # Shared pytest fixtures
    ├── test_titration.py                 # Placeholder test classes (all pass)
    └── fixtures/thesis_data.csv          # 8-row reference dataset (pH, volume, temp)

### Legacy Pascal Source Analysis

File: docs/thesis-reference/5_point_titration_code.pas

The Pascal program (program atct) is approximately 420 lines. It was scanned/OCR'd from a printed thesis, resulting in significant character corruption. Key structural components:

Mathematical functions (lines 41–73):

    log(x) — base-10 logarithm via ln(x)/2.302585093
    tento(y) — 10^y via exp(y*ln(10))
    nue(TDS,dil) — ionic strength: 0.000025 * (TDS/dil - 20)
    logf(mue,ktemp) — Debye-Hückel activity coefficient (heavily OCR-corrupted)
    perH2CO3(ph), perHCO3(ph), perCO3(ph) — carbonate species fractions
    per(ph,pkk) — generic weak acid ionization fraction
    dH2CO3alk(pHf,pHs) — carbonate alkalinity difference between two pH points
    dHAcalk(pHf,pHs) — acetic acid alkalinity difference
    MH2O(vxfi,vxs,pHfi,pHs) — water self-ionization mass balance
    MNH3(pHf,pHs) — ammonia mass balance correction
    MHP04(pHf,pHs) — phosphate mass balance correction

Thermodynamic constants (procedure pK, lines 287–300):

    Temperature-dependent pK calculations for carbonate (pK1, pK2), acetic acid (pKa), ammonia (pKn), phosphate (pKp)
    Activity-corrected versions (pK11, pK22, pKaa, pKnn, pKpp) using logf

Main solver (procedures atctcalculation and atct1, lines 309–360):

    Calculates Ct (total carbonate) from two independent pH pair combinations (pH1-pH2 and pH1-pH4, both using pH3-pH4 as reference)
    Computes Ctcomp = Ct1 - Ct2 as convergence check
    Iteratively adjusts pH by ±0.01 increments (up to 20 iterations) to minimize Ctcomp
    Calculates At (total SCFA) and H2CO3* alkalinity from converged Ct value

UI procedures (lines 77–285): CRT-based text interface — introscreen, display, screen, box, etc. These are NOT converted; replaced by CLI.

Default values (lines 153–169): Reference dataset from the thesis: pH0=7.36, pH1=6.75, pH2=5.95, pH3=5.18, pH4=4.29, Vx1=1.06, Vx2=3.50, Vx3=4.84, Vx4=5.40, ca=0.0728 mol/L, vsundil=10 mL, vsdil=50 mL, temp=21°C, TDS=3300 mg/L.
Existing Python Drafts

docs/thesis-reference/5point_titration.py — Partial translation using math module. Defines log(), ten_to(), perH2CO3(), perHCO3(), perCO3(), and a basic calculate_titration(). Issues: uses global pK1/pK2 variables, hardcodes pK constants instead of computing them from temperature, missing functions (nue, logf, MH2O, MNH3, MHP04, per, dHAcalk), dHAcalk incorrectly calls perHCO3 instead of per, no iterative pH correction.

docs/thesis-reference/5point_titration_ver1.py — Another attempt with input()-based data entry. Similar issues plus incorrect perHCO3 and perCO3 formulas (missing terms in denominators), undefined per() function called by dHAcalk.
Current Python Scaffold

src/titration/core.py — Contains three stub functions:

    find_endpoint() — Uses np.gradient() for derivative-based endpoint detection (not in the Pascal original; appears to be a different algorithm)
    calculate_alkalinity() — Simple placeholder formula, not the full 5-point method
    parse_pascal_equivalent() — Returns a status dict, no actual parsing

src/titration/__init__.py — Exports the three stubs, version 0.1.0.
# Build and Dependencies

pyproject.toml: setuptools-based build, Python ≥ 3.10. Core dependencies: numpy ≥ 1.26.0, scipy ≥ 1.12.0, pandas ≥ 2.1.0. Dev dependencies: pytest ≥ 7.4.0, pytest-cov ≥ 4.1.0, matplotlib ≥ 3.8.0.

requirements.txt: Adds seaborn ≥ 0.13.0 and tqdm ≥ 4.66.0 beyond what pyproject.toml lists.
CI Pipeline

.github/workflows/python-tests.yml: Runs on push/PR to main. Tests across Python 3.10, 3.11, 3.12. Installs via pip install -r requirements.txt then pip install -e ".[dev]". Sets PYTHONPATH=src. Runs pytest tests/ -v (coverage temporarily disabled). Codecov upload step present but effectively no-op since coverage XML is not generated.
# Test Infrastructure

    tests/conftest.py: Defines shared fixtures (sample_ph_data, sample_volume_data, expected_alkalinity) — note these are generic titration curve data, not the specific 5-point method inputs
    tests/test_titration.py: Two test classes with all pass bodies — TestTitrationAlgorithms and TestDataLoading
    tests/fixtures/thesis_data.csv: 8-row CSV with pH/volume/temperature columns — appears to be a general titration curve, not the 5-point specific inputs (which are 5 pH values + 4 volumes + additional parameters)

Key OCR Artifacts in Pascal Source
OCR'd Text	Correct Pascal	Meaning
111Je / Ill.le	mue	Ionic strength (μ)
14OOO	14000	Molecular weight factor for nitrogen
31OOO	31000	Molecular weight factor for phosphorus
stringt50J	string[50]	Pascal string type
tento(−pk11)	tento(-pk11)	OCR en-dash vs minus
Ctc−""	Ctcomp	Ct comparison variable
kt / ktemp	ktemp	Temperature in Kelvin
Wikis 2/2
architecture/architecture-overview.md
Architecture Overview
Target System Architecture

The modernized 5-point titration package is a pure-Python computational library with no external service dependencies. It converts the monolithic Turbo Pascal 4.0 program atct into a modular Python package (titration) with clearly separated concerns, typed data structures, and NumPy-based math.
Package Structure

src/titration/
├── __init__.py       # Public API exports, version
├── models.py         # TitrationInput / TitrationResult dataclasses
├── constants.py      # Thesis default values, conversion factors
├── chemistry.py      # Thermodynamic pK, activity coefficients, species fractions, mass balance
├── solver.py         # Iterative atct solver (atctcalculation + atct1 convergence loop)
├── core.py           # High-level orchestration: run_titration(input) → result
└── cli.py            # argparse CLI entry point

## Module Responsibilities

models.py — Defines the two primary data structures that replace the Pascal program's 30+ global variables:

    TitrationInput: 16 typed fields for pH values, volumes, titrant normality, sample sizes, temperature, TDS, nitrogen, phosphorus
    TitrationResult: computed outputs (H2CO3* alkalinity in mg/L as CaCO3, SCFA in mg/L as acetic acid, systematic pH error, convergence status, intermediate Ct/At values)

constants.py — Conversion factors from the thesis (50000 for CaCO3, 60000 for acetic acid, 14000 for nitrogen, 31000 for phosphorus) and the THESIS_DEFAULTS instance providing the reference dataset from the Pascal default_values procedure.

chemistry.py — All thermodynamic and speciation computations, mapping directly to Pascal functions:

    calculate_ionic_strength(tds, dil) ← Pascal nue(TDS, dil)
    calculate_log_activity(mue, ktemp) ← Pascal logf(mue, ktemp) (Debye-Hückel equation)
    calculate_pk_constants(temperature, tds, dil) ← Pascal pK procedure — returns dict with all raw and activity-corrected pK values
    per_h2co3(ph, pk11, pk22) ← Pascal perH2CO3(ph)
    per_hco3(ph, pk11, pk22) ← Pascal perHCO3(ph)
    per_co3(ph, pk11, pk22) ← Pascal perCO3(ph)
    per(ph, pkk) ← Pascal per(ph, pkk)
    d_h2co3_alk(ph_f, ph_s, pk11, pk22) ← Pascal dH2CO3alk(pHf, pHs)
    d_hac_alk(ph_f, ph_s, pkaa) ← Pascal dHAcalk(pHf, pHs)
    m_h2o(vxfi, vxs, ph_fi, ph_s, vsdil, logf1) ← Pascal MH2O(vxfi, vxs, pHfi, pHs)
    m_nh3(ph_f, ph_s, nt, dil, vsdil, pknn) ← Pascal MNH3(pHf, pHs)
    m_hpo4(ph_f, ph_s, pt, dil, vsdil, pkpp) ← Pascal MHP04(pHf, pHs)

All functions are implemented as pure functions accepting explicit parameters (no global state). They use NumPy scalar functions (np.log10, np.sqrt) so they work with both float and np.ndarray inputs.

solver.py — Implements the iterative convergence algorithm:

    run_solver(inp: TitrationInput) -> TitrationResult ← Pascal atct1 — the main entry point that computes dilution factor, pK constants, then iteratively adjusts all pH values by ±0.01 (up to 20 iterations) until Ct1 ≈ Ct2, then computes At (SCFA) and H2CO3* alkalinity. Works on copies of pH values to avoid mutating the input.
    _atctcalculation() (nested closure) ← Pascal atctcalculation — computes Ct1 and Ct2 from two independent pH pair combinations (pH1-pH2 and pH1-pH4, both referenced against pH3-pH4), returns (m_ct1, ct1, ct2, ct_comp).

core.py — Thin orchestrator with a single entry point run_titration(inp: TitrationInput) -> TitrationResult that delegates to run_solver.

cli.py — Thin argparse wrapper registered as titration console script in pyproject.toml. Accepts input parameters as command-line arguments (with thesis defaults as fallbacks), calls run_titration, formats output.
Key Architectural Patterns
Pure Functions with Explicit Parameters

The Pascal program relies on global mutable state — procedures like perH2CO3(ph) implicitly read global pk11 and pk22. The Python equivalent takes all dependencies as explicit parameters: species_fraction_h2co3(ph, pk11, pk22). This enables independent unit testing and future vectorized batch processing.
Immutable Input, Working Copies for Mutation

The Pascal atct1 procedure mutates global pH values via deltapH during iteration. The Python solver receives an immutable TitrationInput and creates a working copy of pH values to mutate within the iteration loop. The original input is never modified.
Layered Computation

Data flows in one direction through the module layers:

TitrationInput
     │
     ▼
constants.py    (conversion factors)
     │
     ▼
chemistry.py    (pK values, activity coefficients, species fractions)
     │
     ▼
solver.py       (iterative Ct convergence, At calculation, alkalinity)
     │
     ▼
TitrationResult

core.py orchestrates this flow. cli.py is a thin shell over core.py.
NumPy Scalar-Compatible Functions

All mathematical functions use np.log10, np.sqrt, etc. instead of math.log10, math.sqrt. These work identically on scalar float inputs but also accept np.ndarray, meaning batch processing of multiple titrations can be added in the future without rewriting function internals.
Post-Modernization Codebase Structure

5-point-titration/
├── .github/workflows/python-tests.yml
├── pyproject.toml                       # Updated with CLI entry point
├── requirements.txt
├── src/titration/
│   ├── __init__.py                      # Updated exports
│   ├── models.py                        # NEW — dataclasses
│   ├── constants.py                     # NEW — defaults, conversion factors
│   ├── chemistry.py                     # NEW — all thermodynamic/speciation math
│   ├── solver.py                        # NEW — iterative solver
│   ├── core.py                          # REWRITTEN — orchestrator
│   └── cli.py                           # NEW — CLI entry point
├── tests/
│   ├── conftest.py                      # Shared fixtures (thesis_defaults, thesis_pk_constants)
│   ├── test_chemistry.py               # Unit tests for chemistry functions (35 tests)
│   ├── test_solver.py                  # Integration tests for solver (13 tests)
│   ├── test_titration.py              # Parametric and boundary tests (23 tests)
│   └── fixtures/thesis_data.csv         # Reference data
├── .mcode/functional-tests/             # Functional test artifacts from automated testing
└── docs/thesis-reference/               # Original Pascal source and draft translations

# Architecture
## Overview
graph TD
    subgraph Legend
        L1[Modified]:::modified
        L2[Unchanged]:::unchanged
    end

    CLI["cli.py\n(CLI Entry Point)"]:::modified
    CORE["core.py\n(Orchestrator)"]:::unchanged
    SOLVER["solver.py\n(Iterative Solver)"]:::modified
    CHEM["chemistry.py\n(Thermodynamics)"]:::unchanged
    MODELS["models.py\n(Data Structures)"]:::unchanged
    CONST["constants.py\n(Defaults & Factors)"]:::unchanged

    CLI --> CORE
    CORE --> SOLVER
    SOLVER --> CHEM
    SOLVER --> MODELS
    CORE --> MODELS
    CLI --> MODELS
    SOLVER --> CONST
    CONST --> MODELS

    classDef modified fill:#ffd700,stroke:#333,color:#000
    classDef unchanged fill:#f0f0f0,stroke:#333,color:#000


## Technology Choices
Component	Choice	Rationale
Python ≥ 3.10	Language	Dataclass features, union type syntax, match statements
NumPy	Math library	Scalar-compatible functions, future vectorization path
SciPy	Scientific computing	Available for potential optimization extensions
pandas	Data handling	Available for batch I/O and CSV processing
pytest + pytest-cov	Testing	Already configured in CI
argparse (stdlib)	CLI	No additional dependency, sufficient for parameter parsing
setuptools	Build	Already configured in pyproject.toml

All dependencies use permissive open-source licenses (BSD-3-Clause, MIT, PSF).
