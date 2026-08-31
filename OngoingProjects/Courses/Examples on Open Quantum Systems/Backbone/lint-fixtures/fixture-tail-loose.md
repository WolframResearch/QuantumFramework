# Fixture: a tail that passes the loose verification check and must fail the strict one

This fixture reproduces the loose-tail defect: S1's tail says "Verified against the companion battery", which the pre-2026-08 linter accepted (any "verified") but which names none of the typed verification kinds. The linter must raise exactly one ERROR, on S1, plus the forward-application NOTICE for the undeclared companion.

## The setup

**R1 (assumption).** The system couples to a single lossless probe, and nothing else touches it.

## The equation

**S0 (from R1; cited).** The probe turns the setup into a linear input relation between the system variable and the record. *Cited: Alpha and Beta, Journal of Fixtures 1, 1 (2000). Verification: none; the translation is trusted to its source.*

## Consequences

**S1 (from S0).** The record's variance falls inversely with integration time. *Derived in place. Verified against the companion battery, whole node; no shared code.*

## The headline

**H (from S1).** The variance floor is set by the probe alone. Audit: H rests on the contract S0; S0 is realized from R1 by the cited translation. No conjecture anywhere in the closure. *Derived in place from S1. Verified by numerical check, whole node; no shared code.*
