# Fixture: a tail claims a verification but the companion has no anchor for it

This fixture reproduces the manifest-anchor defect: the page declares its companion, S1's tail claims a numerical check, and the companion carries an anchor for H but none for S1. The linter must raise exactly one ERROR, on S1, under rule 5.

<!-- companion: fixture-manifest-companion.md -->

## The setup

**R1 (assumption).** The system couples to a single lossless probe, and nothing else touches it.

## The equation

**S0 (from R1; cited).** The probe turns the setup into a linear input relation between the system variable and the record. *Cited: Alpha and Beta, Journal of Fixtures 1, 1 (2000). Verification: none; the translation is trusted to its source.*

## Consequences

**S1 (from S0).** The record's variance falls inversely with integration time. *Derived in place. Verified by numerical check, whole node, at a single coupling; no shared code.*

## The headline

**H (from S1).** The variance floor is set by the probe alone. Audit: H rests on the contract S0; S0 is realized from R1 by the cited translation. No conjecture anywhere in the closure. *Derived in place from S1. Verified by numerical check, whole node; no shared code.*
