# Lanyon AI: corporate profile (identity, people, funding, history)

Research date: 2026-10-04. Each finding carries an evidence tag. **[REG]** marks a government registry or filing search. **[TECH]** marks a technical record: WHOIS/RDAP, the GitHub API, the Wayback Machine, or timestamps decoded from X post IDs. **[CO]** marks a company self-report (website, blog, press release, company or founder X posts). **[INST]** marks an institutional or academic page. **[DB]** marks a third-party startup database. **[PRESS]** marks third-party press; nearly all Lanyon press is a rewrite of the company's own release. Inferences appear only in the "Inferences" subsections. Unless stated otherwise, every page below was accessed on 2026-10-04.

## 1. Which entity is "Lanyon AI"? (identity, domain, legal name, incorporation, HQ, rule-outs)

### Takeaway
"Lanyon AI" is **Lanyon AI, Inc.** (lanyon.ai), a Princeton, NJ company that calls itself a "fundamental research lab". It came out of stealth on 2026-07-17 and announced a $10.6M round on 2026-08-17. Its AI agent, also called "Lanyon", generates numerical PDE solvers (C code) together with machine-checked correctness proofs (Lean/Rocq). It is the only "Lanyon AI" I found, and it matches "verification of numerical solvers / scientific simulation software". New Jersey lists it as a *foreign* (out-of-state) for-profit corporation, registered in NJ on 2026-08-14. Its home state and original incorporation date could not be established.

### Cited Findings
**Legal name and registration**
- The website footer reads "© 2026 Lanyon AI, Inc." [CO], 2026-10-04 — [lanyon.ai](https://lanyon.ai/); [Company page](https://lanyon.ai/company/)
- The press release (2026-08-17) closes with "SOURCE Lanyon AI Inc." and a "Lanyon AI Inc." signature block [CO] — [PR Newswire, 2026-08-17](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- **NJ registry [REG]**, from a NJ Division of Revenue business-name search for "Lanyon" run 2026-10-04:
  - Name: "LANYON AI INC."
  - Entity ID 0451511506
  - City: PENNINGTON
  - Type: FR, which the site expands as "Foreign For-Profit Corporation"
  - Date: 8/14/2026, shown in a column headed "Incorporated Date"
  - Source: [NJ Business Name Search](https://www.njportal.com/DOR/BusinessNameSearch/Search/BusinessName)
- The American Physical Society's donor list for the William D. Dorland Prize endowment names "Lanyon AI, Inc." among organizational donors. This is independent third-party confirmation of the legal name [INST], accessed 2026-10-04 — [APS Dorland Prize page](https://www.aps.org/about/support/support-dorland-prize-computational-plasma)
- **SEC EDGAR [REG]**, checked 2026-10-04:
  - No registrant is named "Lanyon AI".
  - A full-text search for "Lanyon AI" over 2026-01-01 to 2026-10-04 returned 0 hits.
  - The only EDGAR company beginning "Lanyon" is the unrelated **Lanyon, Inc.** (CIK 0001476844), of Irving, TX, incorporated in DE, which filed a Form D on 2009-11-24.
  - Sources: [EDGAR company search "lanyon"](https://www.sec.gov/cgi-bin/browse-edgar?company=lanyon&owner=exclude&action=getcompany); [EDGAR full-text search](https://efts.sec.gov/LATEST/search-index?q=%22Lanyon%20AI%22&dateRange=custom&startdt=2026-01-01&enddt=2026-10-04)
- A UK Companies House search for "lanyon ai" (2026-10-04) returned no Lanyon AI entity. The results are unrelated: Lanyon Bowdler LLP, Belfast "Lanyon Place/Quay" property companies, and similar [REG] — [Companies House search](https://find-and-update.company-information.service.gov.uk/search/companies?q=lanyon+ai)

**Headquarters**
- The press release (2026-08-17) gives "100 Overlook Center, Suite 2145, Princeton, NJ 08540" [CO] — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- 100 Overlook Center (2nd floor), Princeton 08540, is a **Regus** serviced-office/coworking centre. Regus sells private offices, coworking, "virtual office plans" and a "business address" product there [third-party listing], accessed 2026-10-04 — [Regus Overlook Center](https://www.regus.com/en/us/433)
- Dealroom gives HQ Princeton, US, and founding year 2026 [DB], accessed 2026-10-04 — [Dealroom](https://dealroom.co/companies/lanyon-ai/)

**Domain and online footprint**
- **lanyon.ai domain [TECH]**, per RDAP:
  - Registered 2026-06-05T02:48:42Z
  - Registrar: InternetX
  - Expires 2028-06-05
  - Last changed 2026-08-12
  - Registrant masked by "PrivateName Services Inc." (Vancouver, BC)
  - Source: [RDAP record](https://rdap.identitydigital.services/rdap/domain/lanyon.ai)
- **GitHub organization `lanyonai` [TECH]**:
  - Display name "Lanyon AI Inc"
  - Location "United States of America"
  - Description "Formal Verification for a Computable Universe"
  - Created 2026-04-04T17:24:33Z
  - Source: [GitHub org](https://github.com/lanyonai); [GitHub API](https://api.github.com/orgs/lanyonai)
- **X account @lanyon_ai [TECH]** ("Lanyon AI"):
  - User ID 2064912680138199040, which decodes as a Twitter snowflake to account creation on about 2026-06-11.
  - First post: 2026-07-17 13:04 UTC.
  - Source: [@lanyon_ai](https://x.com/lanyon_ai)

**What the company does (only as far as needed to identify it)**
- Self-description: "Lanyon AI, a fundamental research lab developing a new kind of scientific and technical AI" [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- Homepage tagline: "Formal Verification for a Computable Universe … building the formally verified substrate connecting AI to the physical world." Its "Specification → Implementation + Proof" pitch says a "neurosymbolic compiler generates optimized implementations and proofs, deterministically, from the same DSL source" [CO], 2026-10-04 — [lanyon.ai](https://lanyon.ai/)
- The public repos are each described as "End-to-end formally verified solvers for…" a named equation system: advection-diffusion, Maxwell, electrostatic Vlasov, Burgers, general-relativistic Maxwell, compressible Euler, ideal MHD [TECH/CO], 2026-10-04 — [github.com/lanyonai](https://github.com/lanyonai)

**No second "Lanyon AI"**
- Searches for "Lanyon AI", "LANYON AI INC" and Lanyon + AI variants on 2026-10-04 returned only this company and its press syndications. No second "Lanyon AI" was found.

**Unrelated "Lanyon" entities, ruled out**
- **Lanyon Solutions, Inc.** made SaaS for meetings, events and travel and "was merged into Cvent in 2017" — [Wikipedia: Lanyon Solutions](https://en.wikipedia.org/wiki/Lanyon_Solutions)
  - After acquiring Cvent, Vista Equity Partners "would merge Cvent with Lanyon, another meetings-technology firm owned by Vista" — [Wikipedia: Cvent](https://en.wikipedia.org/wiki/Cvent); [Cvent press release "Cvent and Lanyon Announce Merger"](https://www.cvent.com/en/press-release/cvent-and-lanyon-announce-merger)
  - Its registry traces are distinct from Lanyon AI's [REG]:
    - NJ lists "LANYON SOLUTIONS, INC.", type FR, Tysons Corner, 8/31/2000 — [NJ search](https://www.njportal.com/DOR/BusinessNameSearch/Search/BusinessName)
    - EDGAR lists Lanyon, Inc. (Irving TX, DE), with a Form D in 2009 — [EDGAR](https://www.sec.gov/cgi-bin/browse-edgar?action=getcompany&CIK=0001476844)
    - Companies House lists the overseas company "LANYON, INC." FC031137 — [Companies House](https://find-and-update.company-information.service.gov.uk/company/FC031137)
- **Lanyon Bowdler** is a full-service law firm in Shropshire and Herefordshire, SRA number 534828 — [lblaw.co.uk](https://www.lblaw.co.uk/); [SRA register](https://www.sra.org.uk/consumers/register/organisation/?sraNumber=534828)
  - Companies House lists LANYON BOWDLER LLP OC351948, incorporated 31 Jan 2010, Shrewsbury [REG] — [Companies House search](https://find-and-update.company-information.service.gov.uk/search/companies?q=lanyon+ai)
- **Lanyon Jekyll theme**: `poole/lanyon`, "A content-first, sliding sidebar theme for Jekyll", repo created 2013-12-28 [TECH] — [GitHub poole/lanyon](https://github.com/poole/lanyon)
- **"Lanyon" Markdown web server in Go** (HN, 2014), named after Dr. Hastie Lanyon in *Jekyll and Hyde* — [HN item](https://news.ycombinator.com/item?id=7713941)
- **Places and people**:
  - Lanyon Quoit, a dolmen in Cornwall — [Wikipedia](https://en.wikipedia.org/wiki/Lanyon_Quoit)
  - Sir Charles Lanyon, a 19th-century Belfast architect; Belfast's "Lanyon Place" and "Lanyon Quay" companies appear in Companies House results — [Wikipedia](https://en.wikipedia.org/wiki/Charles_Lanyon)
- **Other NJ "Lanyon" entities**, all unrelated: Lanyon & Irvin LLC (Ocean Grove, 2010), Lanyon Management LLC (Wayside, 2005) and The Lanyon Group, Inc. (Weehawken, 2012) [REG] — [NJ search](https://www.njportal.com/DOR/BusinessNameSearch/Search/BusinessName)

### Inferences
- **Home state.** NJ classifies the company as a *foreign* corporation, so it was incorporated outside New Jersey. Delaware is the most likely home state, since it is the default for VC-backed US startups, but this is unconfirmed.
- **When it was incorporated.** Probably in the first half of 2026. Three facts point there:
  - The GitHub org named "Lanyon AI Inc" was created on 2026-04-04, though its display name may have been set later; the org was last updated 2026-07-13.
  - Gorard writes that angel Siqi Chen was among the first people he spoke to "after we incorporated the company", before the round closed in July 2026.
  - The NJ foreign registration (2026-08-14) came three days before the press release, which looks like housekeeping ahead of the announcement.
- **Two addresses.** The NJ record shows Pennington, a town near Princeton, while the press release shows a Princeton Regus suite. The NJ filing probably uses a founder's or agent's address. The Regus suite points to a small serviced office, not a purpose-built lab.

### Gaps
- **State and date of incorporation.**
  - Delaware's ICIS entity search is CAPTCHA-gated and could not be scripted.
  - The OpenCorporates API needs a token, and its web search returned no parseable results.
  - No Form D exists to give either fact.
  - The NJ portal shows only the summary row; the full NJ status report, with registered agent and home jurisdiction, is a paid product.
- The PitchBook profile [ID 1479893-41](https://pitchbook.com/profiles/company/1479893-41) returned HTTP 403, so its founding date and legal-entity fields were not seen.

## 2. Founders and key team

### Takeaway
All three co-founders came from the Princeton Plasma Physics Laboratory (PPPL) Gkeyll computational-plasma group:
- **Jonathan Gorard (CEO)**: applied mathematician and co-founder of the Wolfram Physics Project. He holds a Cambridge MPhil in Scientific Computing, has a background in Wolfram Language automated theorem proving, and joined Princeton as a Research Software Engineer in Feb 2024.
- **Ammar Hakim (CTO)**: PhD in aeronautics/astronautics, University of Washington, 2006; previously at Tech-X. At PPPL he was Principal Research Physicist, Deputy Head of Computational Sciences and lead of Gkeyll.
- **James "Jimmy" Juno (Chief Scientist)**: BS Rice 2014; PhD UMD 2020 under Bill Dorland; PPPL staff research physicist from 2022.

Their expertise covers:
- numerical analysis: discontinuous Galerkin and finite-volume methods for hyperbolic and kinetic PDEs
- HPC and GPUs
- numerical relativity
- formal verification and automated theorem proving, which is mainly Gorard's

Gorard's quantum work (ZX-calculus, Wolfram Language quantum functionality) predates Lanyon. No advisors, board members or non-founder employees are named publicly.

### Cited Findings
**Company-stated roles** [CO], 2026-10-04 — [Company page](https://lanyon.ai/company/):
- **Gorard**, "Co-founder & CEO": "Applied mathematician: responsible for orchestrating Lanyon's overarching R&D mission, and for the design of all formal, symbolic, and mathematical aspects."
- **Hakim**, "Co-founder & CTO": "Algorithm alchemist: responsible for the design and implementation of all of Lanyon's key numerical and computational algorithms."
- **Juno**, "Co-founder & Chief Scientist": "Computational physicist: responsible for closing the gap between mathematical and algorithmic implementation and production scientific problems."
- The company says its founders bring "many decades of pioneering research … spanning automatic code-generation, automated theorem-proving, scientific AI, and the computational foundations of physics."

**Press-release descriptions** [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html):
- The three are "all formerly of Princeton University and/or the Princeton Plasma Physics Laboratory".
- Gorard is "an award-winning applied mathematician, known previously for co-founding the Wolfram Physics Project with Stephen Wolfram".
- Hakim is "a world-leading computational physicist, with deep expertise in fluid mechanics, nuclear fusion, and aerospace engineering".
- Juno is "a leading plasma physicist".
- The team has "over five decades of expertise between them".

**Jonathan Gorard (CEO)**
- **Princeton Research Computing bio** [INST], Wayback capture 2026-02-06 — [Princeton RC (archived)](https://web.archive.org/web/20260206171050/https://researchcomputing.princeton.edu/about/people-directory/jonathan-gorard):
  - Position: "Research Software Engineer II", PPPL
  - "Background: MPhil In Scientific Computing (University of Cambridge)"
  - "joined the RSE group at Princeton in February 2024 … part of the Gkeyll team, led by Ammar Hakim", working on numerical general relativity, adaptive mesh refinement and curvilinear meshing
  - Previously "a graduate student at the University of Cambridge" who "held various positions (including as a research fellow and director of research) at Cardiff University, Wolfram Research and the Wolfram Institute"
  - Interests: numerical relativity and applied category theory
- **Wolfram Physics Project page** [INST], undated, pre-2024 framing — [wolframphysics.org](https://www.wolframphysics.org/pages/people/jonathan-gorard/):
  - Title "Associate Director of Research & Academics"
  - "consultant mathematician for Wolfram Research (leading the development of the Wolfram Language's automated theorem proving, axiomatic mathematics, quantum computing and discrete-state quantum mechanics functionality)"
  - Co-founder of the Wolfram Physics Project
  - Location UK; affiliation "University of Cambridge and Wolfram Research"
  - The Wolfram Institute page carries the same bio — [wolframinstitute.org](https://wolframinstitute.org/people/jonathan-gorard)
- **Wolfram Summer Research Institute 2017** [INST] — [Wolfram alumni page](https://education.wolfram.com/summer-research-institute/alumni/2017/gorard/):
  - Described as "a student and junior researcher in the department of mathematics at King's College London"
  - Project: "Automated Theorem Proving for Equational Logic", extending FullSimplify to generate proofs
- **Google Scholar** (2026-10-04): 443 citations, h-index 11, i10-index 12; affiliation still shown as Princeton University — [Scholar](https://scholar.google.com/citations?user=ItG_Nz0AAAAJ)
- **arXiv record** [TECH], 23 papers, 2016–2026 — [arXiv API query](https://export.arxiv.org/api/query?search_query=au:Gorard&sortBy=submittedDate&sortOrder=descending&max_results=45):
  - Earliest: "Uniqueness Trees: A Possible Polynomial Approach to the Graph Isomorphism Problem" (2016)
  - Wolfram-model physics, 2020–2023
  - ZX-calculus: "ZX-Calculus and Extended Wolfram Model Systems II: Fast Diagrammatic Reasoning with an Application to Quantum Circuit Simplification" (2021) — [arXiv:2103.15820](https://arxiv.org/abs/2103.15820)
  - "Computational General Relativity in the Wolfram Language using Gravitas I/II" (2023/2024) — [arXiv:2308.07508](https://arxiv.org/abs/2308.07508)
  - "Applied Category Theory … using Categorica I" (2024) — [arXiv:2403.16269](https://arxiv.org/abs/2403.16269)
  - "Quantum Cellular Automata, Black Hole Thermodynamics, and the Laws of Quantum Complexity" (2019) — [arXiv:1910.00578](https://arxiv.org/abs/1910.00578)
- **Precursors to Lanyon's technology** [TECH/CO]:
  - **Shock with Confidence.** "Shock with Confidence: Formal Proofs of Correctness for Hyperbolic Partial Differential Equation Solvers" (Gorard & Hakim), arXiv, 2025-03-18 — [arXiv:2503.13877](https://arxiv.org/abs/2503.13877)
    - The abstract describes "a new formal verification pipeline for such algorithms in Racket". It builds a bespoke hyperbolic PDE solver, generates "low-level C code which verifiably implements that solver", and produces formal proofs.
    - Gorard's X thread (2025-03-19) called it "the first automated theorem-proving framework for (hyperbolic) PDE solvers" (1,750 likes at retrieval) — [X post](https://x.com/getjonwithit/status/1902158541839856071)
  - **BEACONS.** "BEACONS: Bounded-Error, Algebraically-Composable Neural Solvers for PDEs" (Gorard, Hakim, Juno), 2026-02-16. It covers "formally-verified neural network solvers for PDEs, with rigorous convergence, stability, and conservation properties" — [arXiv:2602.14853](https://arxiv.org/abs/2602.14853)
  - **Joint numerical-relativity papers:**
    - "A Tetrad-First Approach to Robust Numerical Algorithms in General Relativity" (Gorard, Hakim, Juno, TenBarge, 2024) — [arXiv:2410.02549](https://arxiv.org/abs/2410.02549)
    - "Beyond GRMHD" (Gorard, Juno, Hakim, 2025) — [arXiv:2510.26019](https://arxiv.org/abs/2510.26019)
  - **Inverse problems.** "Improved Dimensionality Reduction for Inverse Problems in Nuclear Fusion and High-Energy Astrophysics" (Gorard, Hakim, Hong Qin, Kyle Parfrey, Shantenu Jha, 2025) — [arXiv:2505.03849](https://arxiv.org/abs/2505.03849)
  - None of the three founders has an arXiv paper dated after BEACONS (2026-02-16), and none carries a Lanyon affiliation, as of 2026-10-04 — [arXiv API query](https://export.arxiv.org/api/query?search_query=au:Gorard&sortBy=submittedDate&sortOrder=descending&max_results=45)

**Ammar Hakim (CTO)**
- **PPPL bio** [INST], Wayback capture 2026-08-02 (still live then) — [PPPL (archived)](https://web.archive.org/web/20260802195838/https://www.pppl.gov/people/ammar-hakim):
  - Title: "Principal Research Physicist"
  - "principal computational physicist and Deputy Head of the Computational Sciences Department … where he leads the Applied and Computational Mathematics group"
  - "Ph.D. in aerospace engineering from the University of Washington … high-resolution numerical methods for two-fluid plasma simulations"
  - Previously at "Tech-X Corporation, where he led projects in computational electromagnetism and plasma fluid modeling"
- **Dissertation** [INST]: "High Resolution Wave Propagation Schemes for Two-Fluid Plasma Simulations", University of Washington, 2006; supervisory committee chair Uri Shumlak — [UW thesis PDF](https://www.aa.washington.edu/sites/aa/files/research/cpdlab/docs/PhDthesis_hakim.pdf)
- **GitHub bio** (2026-10-04): "CTO, Lanyon AI (lanyon.ai) … Algorithm Alchemist and Lead Physicist for the Gkeyll Group"; company field "Lanyon AI" [TECH/self] — [github.com/ammarhakim](https://github.com/ammarhakim)
  - His LinkedIn profile is indexed with the headline "Ammar Hakim - Lanyon AI" (page not retrievable) — [LinkedIn](https://www.linkedin.com/in/ammar-hakim-6785b825/)
- **Selected HPC/numerics papers** [self-CV] — [Juno CV](https://spacalum.rice.edu/bios/james_juno_cv.pdf):
  - Hakim & Juno, "Alias-Free, Matrix-Free, and Quadrature-Free Discontinuous Galerkin Algorithms for (Plasma) Kinetic Equations", SC20 (IEEE, 2020)
  - Hakim, Francisquez, Juno, Hammett, "Conservative discontinuous Galerkin schemes for nonlinear Dougherty–Fokker–Planck collision operators", J. Plasma Phys. 86(4), 2020
- **PPPL roles, per Lanyon's blog** [CO], Sep 2026 — [Lanyon blog: Dorland Prize](https://lanyon.ai/blog/dorland/):
  - Bill Dorland "served as the department head of the Princeton Plasma Physics Laboratory's (PPPL) Computational Sciences Department while Ammar was the deputy head and Jimmy was a staff research physicist."

**James "Jimmy" Juno (Chief Scientist)**
- **CV** [self] — [Juno CV, Rice SPAC alumni](https://spacalum.rice.edu/bios/james_juno_cv.pdf):
  - Education:
    - "2014 BS in Computational Physics, Rice University"
    - "2020 PhD in Physics, University of Maryland, College Park", thesis "A Deep Dive into the Distribution Function: Understanding Phase Space Dynamics using Continuum Vlasov–Maxwell simulations", advisers William Dorland and Jason TenBarge
  - Employment:
    - Assistant Research Scientist, PPPL, Jun–Sep 2020
    - NSF AGS Postdoctoral Fellow, University of Iowa, Sep 2020–Jan 2022
    - Staff Research Physicist, PPPL, Jan 2022–
  - Awards:
    - NASA Earth and Space Science Fellowship, 2017 ($135,000)
    - NSF AGS Postdoctoral Fellowship, 2020 ($190,000)
    - DOE National Undergraduate Fellowship in Plasma Physics, 2013
- **PPPL bio** [INST], Wayback capture 2026-01-06 — [PPPL (archived)](https://web.archive.org/web/20260106231441/https://www.pppl.gov/people/james-jimmy-juno):
  - "Staff Research Physicist in the Computational Sciences Department … core development team for the Gkeyll simulation framework", maintaining "the one-of-a-kind multi-species, continuum Vlasov-Maxwell solver"
  - Was developing relativistic-plasma capability for "PPPL's extreme astrophysics initiative"
  - Postdoc mentor: Greg Howes, Iowa
- **Google Scholar** (2026-10-04): 1,083 citations, h-index 16, i10-index 30; affiliation still PPPL — [Scholar](https://scholar.google.com/citations?user=5xPBjHkAAAAJ)
- **Thesis on arXiv** — [arXiv:2005.13539](https://arxiv.org/abs/2005.13539)
- **Recent paper:** "Modeling of Relativistic Plasmas with a Conservative Discontinuous Galerkin Method" (2026-02-19) — [arXiv:2602.17487](https://arxiv.org/abs/2602.17487)
- **His Dorland blog post** (Sep 2026) adds [CO] — [Lanyon blog: Dorland Prize](https://lanyon.ai/blog/dorland/):
  - He arrived at UMD in summer 2014 and was earlier "an intern at PPPL (with Ammar!)".
  - He defended his PhD by Zoom during COVID.
  - He then moved to "a postdoc and later a staff research physicist at PPPL".
  - He names collaborator Sasha Philippov.
- **GitHub** (`JunoRavin`) still lists "Princeton Plasma Physics Laboratory" as company, 2026-10-04 [TECH] — [github.com/JunoRavin](https://github.com/JunoRavin)

**Wider team, advisors and press contact**
- The company page lists only the three founders. Its hiring blurbs are generic (see §4) [CO], 2026-10-04 — [Company page](https://lanyon.ai/company/)
- No advisors or board members are named on the site, in the press release, or in the accessible databases (Dealroom, Equilar).
- Equilar lists Gorard as "Co-founder & Chief Executive Officer" with a start date of 08/20/2026. That is almost certainly the date Equilar captured the data, not his real start [DB] — [Equilar](https://people.equilar.com/bio/person/jonathan-gorard-lanyon-ai/80919672)
- The press-release contact is "Brian Golden", reachable at contact@lanyon.ai (decoded from the page's obfuscated mailto) [CO], 2026-08-17. No public profile tying a Brian Golden to Lanyon was found — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)

### Inferences
- **Who covers what.** Gorard supplies formal methods, symbolic computation and automated theorem proving (Wolfram, 2017 onward). Hakim and Juno supply numerical analysis of hyperbolic and kinetic PDEs (discontinuous Galerkin, finite volume), plasma physics, and HPC/GPU code; Gkeyll is their production codebase. The fit to "numerical solver verification" is direct.
- **Quantum.** Quantum computing ties appear only in Gorard's Wolfram-era work (ZX-calculus papers; leading Wolfram Language quantum-computing functionality). No Lanyon material lists quantum among its targets.
- **Origin.** In substance, Lanyon looks like a spin-out of the Gkeyll group's 2025–26 verification research (Shock with Confidence, BEACONS), now commercialized. No public document describes any formal licence or IP arrangement with Princeton or PPPL.
- **"Over five decades" of expertise.** This is roughly consistent with the record only if counted generously: Hakim since his early-2000s doctoral work (PhD 2006), Juno since undergraduate research around 2012, and Gorard publishing since 2016.
- **"Decades of formal methods expertise".** Dimension's quote overstates the formal-methods depth. In the public record, that expertise is mainly Gorard's (Wolfram theorem proving from 2017; Shock with Confidence, 2025).

### Gaps
- **Departure dates** from Princeton and PPPL are unconfirmed.
  - PPPL still listed Hakim as "Principal Research Physicist" on 2026-08-02.
  - Juno's GitHub and Scholar still show PPPL.
  - Gorard's Scholar still shows Princeton.
- **Gorard's highest degree.** His Princeton bio lists only the Cambridge MPhil, and I found no source for a completed doctorate.
- **Gorard's award.** No award was found to support "award-winning".
- **Hakim's Scholar metrics.** The Scholar author search did not parse.
- **LinkedIn details** for all founders could not be retrieved (blocked).
- **Advisors, board composition, and any non-founder staff** remain unknown.

## 3. Funding (rounds, investors, accelerators, grants, valuation)

### Takeaway
One round is disclosed: a **$10.6M "initial fundraising round"**. **Dimension** led it, with **Industrious Ventures** and angel **Siqi Chen**. Per the company it closed in **July 2026**; it was announced on **2026-08-17**. The company has disclosed neither the instrument (SAFE or priced) nor a valuation. Dealroom alone shows a $42M post-money valuation, which is unverified. I found no accelerator, grant or government contract, and no SEC Form D on EDGAR as of 2026-10-04.

### Cited Findings
**The round, per the company**
- Press release [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html); syndicated by [Yahoo Finance](https://finance.yahoo.com/technology/ai/articles/lanyon-ai-emerges-stealth-build-070000418.html) and [AOL](https://www.aol.com/articles/lanyon-ai-emerges-stealth-build-070000000.html):
  - The company "has emerged from stealth following a $10.6 million initial fundraising round led by Dimension, with participation from Industrious Ventures."
  - The subhead says the round "backs a team of world-leading Princeton mathematicians and physicists".
- Blog "Introducing our Partners" (Gorard, Aug 2026) [CO] — [Lanyon blog: fundraising](https://lanyon.ai/blog/fundraising/):
  - **Closing and investors.** "Last month, right around the time we emerged from stealth, Lanyon AI officially closed its initial fundraising round, led by Dimension, and with participation from Industrious Ventures and angel investor Siqi Chen."
  - **Dimension contacts.** It names Simon Barnett and Zavain Dar of Dimension. Dimension's "initial questions were not about run rates or market share ('those things will come,' they said), but about research direction."
  - **Industrious contacts.** It names Taylor Sargent of Industrious, who offered "contact with physical reality", meaning introductions in "aerospace, nuclear, and defense". Industrious's Alexandra Johnson "helped us coordinate our first press release".
  - **First investor.** "Siqi Chen was one of the very first people I spoke to after we incorporated the company, and also the first person to invest in us." His money let the company "lease our first office space, buy our first computers, and hire the professional legal counsel that ultimately guided us through the remainder of our fundraising process."
- The press release quotes "Simon Barnett, Partner and Head of Research at Dimension". He frames the bar as "flight controls, nuclear systems, or simulating chip tape-outs" [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)

**The investors**
- **Dimension** [PRESS], 2026-07-21 — [TechCrunch](https://techcrunch.com/2026/07/21/dimension-capitals-800m-third-fund-shows-the-intersection-of-science-and-compute-is-booming/); [Biopharma Dive (Fund III)](https://www.biopharmadive.com/news/dimension-biotech-startups-third-fund-tech-ai/825748/); [Biopharma Dive (launch)](https://www.biopharmadive.com/news/dimension-biotech-venture-firm-tech-startups/640938/):
  - A science-and-compute VC firm founded about 2022 by Zavain Dar and Adam Goulburn (ex-Lux Capital) and Nan Li (ex-Obvious Ventures).
  - It announced an **$800M Fund III on 2026-07-21**, 60% larger than its $500M Fund II.
  - It has about $1.65B under management and 35 portfolio companies.
- **Industrious Ventures** [DB] — [Superscout profile](https://superscout.co/investor/industrious-ventures); [Industrious team](https://industrious.vc/team/):
  - A deep-tech VC founded in 2019, based in Denver and Austin.
  - It invests in aerospace, energy, national security, manufacturing and related sectors.
- **Taylor Sargent** is an Industrious partner focused on aerospace and national security. He was previously a Booz Allen lead scientist supporting DARPA and NASA Goddard [DB] — [NFX Signal](https://signal.nfx.com/investors/taylor-sargent)
- **Siqi Chen** is co-founder and CEO of Runway, a finance-planning software company, and a prolific angel investor [PRESS/DB] — [Cognitive Revolution podcast](https://www.cognitiverevolution.ai/building-an-intelligent-business-os-with-runway-ceo-siqi-chen/); Lanyon's blog links his [Mercury investor-database profile](https://mercury.com/investor-database/siqi-chen)

**Third-party databases**
- **Dealroom** [DB], 2026-10-04 — [Dealroom](https://dealroom.co/companies/lanyon-ai/):
  - "August 2026: $11M raised (Early VC round)"
  - "Valuation: $42M (post-funding)"
  - Investors: Dimension (lead) and Industrious Ventures
- **PitchBook** has a profile (1479893-41) but returned 403 [DB] — [PitchBook](https://pitchbook.com/profiles/company/1479893-41)
- **Gaebler.com** lists Lanyon AI as a "Funded Company"; the page timed out [DB] — [Gaebler](https://www.gaebler.com/Funded-Company-399FADDA-C768-4188-A602-4920E4E40442-Lanyon-AI)

**Filings, accelerators, grants and philanthropy**
- **Form D.** An EDGAR full-text search for "Lanyon AI" (2026-01-01 to 2026-10-04) found nothing, so no Form D has been filed under that name [REG] — [EDGAR FTS](https://efts.sec.gov/LATEST/search-index?q=%22Lanyon%20AI%22&dateRange=custom&startdt=2026-01-01&enddt=2026-10-04)
- **Accelerators.** No Y Combinator, Entrepreneur First or Wellfound listing was found, and no company material mentions an accelerator — search of [ycombinator.com / joinef.com / wellfound.com](https://www.ycombinator.com/companies), 2026-10-04.
- **Grants.** No company material mentions a grant or government contract. The SBIR.gov public API returned 403, so SBIR/STTR could not be checked directly. UKRI and CORDIS are not relevant to a US entity with no UK or EU presence.
- **Philanthropy (an outflow).** Lanyon "contributed to the permanent endowment of the William D. Dorland Prize in Computational Plasma Physics" (amount undisclosed) [CO], Sep 2026 — [Lanyon blog](https://lanyon.ai/blog/dorland/)
  - APS lists "Lanyon AI, Inc." among 47 donors to a $150,000 campaign that stood at $130,518.79 when accessed [INST], 2026-10-04 — [APS](https://www.aps.org/about/support/support-dorland-prize-computational-plasma)

### Inferences
- A $10.6M first round for a three-founder research lab is consistent with a large seed. "Initial" suggests the first institutional round. The instrument is unknown.
- Dealroom's $42M post-money valuation and "$11M" (a rounding of $10.6M) may be Dealroom estimates. Treat them as unverified.
- The round closed in the same month Dimension announced its Fund III. It cannot be determined which Dimension fund invested.
- Industrious is positioned as the go-to-market partner into aerospace, defense and nuclear, consistent with its thesis and Sargent's DARPA background. Dimension is the research-oriented lead.
- **On the missing Form D.** A Form D is normally filed within 15 days of the first sale in a Regulation D offering. If the round closed in July 2026, it would have been due around August, so its absence by 2026-10-04 is notable. It is not conclusive: some issuers rely on §4(a)(2) without filing, some file late, and EDGAR indexing can lag.

### Gaps
- Instrument (SAFE or priced equity), price and valuation; investor-by-investor amounts; other angels; board seats; pro-rata or option-pool details.
- Any grants or government contracts (SBIR, DOE, DARPA): none found and the SBIR API was inaccessible.
- PitchBook and Crunchbase deal records (not accessible).

## 4. Size and traction signals (headcount, hiring, customers, press, social)

### Takeaway
As of 2026-10-04 the company is at a very early stage:
- three named people (the founders) and no named employees, customers or design partners
- a Regus serviced-office suite
- generic open roles for research scientists and engineers
- public output of 9 GitHub repos, 9 research notes (including self-published benchmarks) and 4 blog posts
- press limited to one PR Newswire release and its rewrites
- moderate social reach

### Cited Findings
**People and hiring**
- The company page names only the three founders [CO], 2026-10-04 — [Company page](https://lanyon.ai/company/)
- Open roles, both via a "Come and join us" mailto to contact@lanyon.ai:
  - **"Research Scientists"**: "outstanding researchers in the fields of applied mathematics, formal verification, computational physics, and artificial intelligence"
  - **"Research Engineers"**: "exceptional software engineers and applied researchers, to help us build the next generation of research tooling and AI infrastructure"
- No Lanyon AI LinkedIn company page was found at the obvious slug (linkedin.com/company/lanyon-ai returned 404) [TECH]. No postings were found on job boards (searches, 2026-10-04).
- Dealroom gives "<100" employees, which says nothing useful [DB] — [Dealroom](https://dealroom.co/companies/lanyon-ai/)

**Office**
- A Regus serviced-office suite at 100 Overlook Center, Suite 2145 [CO/third-party] — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html); [Regus](https://www.regus.com/en/us/433)
- The angel money paid for "our first office space" and "our first computers" [CO] — [Lanyon blog: fundraising](https://lanyon.ai/blog/fundraising/)

**GitHub** [TECH], 2026-10-04 — [github.com/lanyonai](https://github.com/lanyonai)
- The org has 139 followers, 0 public members and 9 public repos.
- Repos with stars and creation dates:

| Repo | Stars | Created | Notes |
|---|---|---|---|
| AdvectionDiffusion | 73 | 2026-07-17 | updated 2026-09-09 |
| MaxwellEquations | 44 | 2026-07-20 | |
| GeneralRelativisticMaxwell | 20 | 2026-07-27 | |
| BurgersEquation | 14 | 2026-07-23 | updated 2026-09-09 |
| ElectrostaticVlasov | 12 | 2026-07-21 | |
| CompressibleEuler | 12 | 2026-07-29 | |
| IdealMHD | 9 | 2026-08-27 | |
| gkeyll | 1 | 2026-08-24 | fork of `gkeyllorg/gkeyll` (MIT); last push 2026-09-23 |
| LanyonScripts | 1 | 2026-08-27 | |

- Licences (GitHub API, 2026-10-04): 7 of the 8 Lanyon-authored repos are MIT-licensed, as is the gkeyll fork. AdvectionDiffusion, the first and most-starred repo, has no licence file.
- The published artifacts are the *outputs*: generated C code, Lean proofs and screen captures. The generator is described as "our own proprietary domain-specific language (DSL)" and "our proprietary AI agent" [CO] — [Lanyon blog](https://lanyon.ai/blog/welcome/)
- Commits on the sampled solver repos come almost entirely from `JonathanGorard`, plus one merge by Juno and one small external fix.

**Self-reported output**
- The IdealMHD README reports "~434 seconds for Lanyon to generate everything", "49,317 lines of Lean 4 code", "31,964 lines of formally verified C code", and "270 total definitions and 156 total theorems" [CO], 2026-10-04 — [IdealMHD README](https://github.com/lanyonai/IdealMHD)
- **Research notes** [CO], 2026-10-04 — [Research notes](https://lanyon.ai/research/):
  - 9 published notes plus one "in preparation" (the Godunov-Peshkov-Romenskii continuum-mechanics model).
  - Two are titled "Benchmarking Lanyon against Frontier Models" (linear and nonlinear PDEs).
  - None is dated or attributed to an author.

**Customers**
- No customers, design partners or pilots are named anywhere. The press release names only target industries: "aerospace engineering, space and atmospheric propulsion, and nuclear energy" [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- The founders acknowledge that persuading "aerospace, nuclear, and defense" buyers "was going to be a tall order" for "three nerdy theorists fresh out of academia" [CO], Aug 2026 — [Lanyon blog: fundraising](https://lanyon.ai/blog/fundraising/)

**Press coverage** [PRESS]
- The PR Newswire release (2026-08-17) was syndicated by Yahoo Finance and AOL. Rewrites appeared on the same day or shortly after:
  - [AI Journal](https://aijourn.com/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing/)
  - [TheSaaSNews ("Lanyon AI Raises $10.6M in Funding")](https://www.thesaasnews.com/news/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing/)
  - [TechEdge AI](https://techedgeai.com/lanyon-ai-emerges-with-10-6m-bet-on-provably-correct-scientific-ai/)
  - [Aerospace Trends](https://www.aerospace-trends.com/lanyon-ai-announces-stealth-exit-and-advances-in-high-precision-scientific-computing/), dated 2026-08-17, with no facts beyond the release
- Searches on 2026-10-04 found no TechCrunch, Axios, Forbes, Business Insider or Sifted coverage of Lanyon, and no Hacker News thread.
- A Digg item, "Jonathan Gorard releases Lanyon, an AI agent that uses …", now returns 404 [aggregator] — [Digg](https://digg.com/tech/31wa8sqc)
  - Search snippets attribute to it the claim that Lanyon "can one-shot a provably-correct C implementation of a complex 3D nonlinear PDE solver (e.g. resistive MHD), and synthesize a ~10k-line Lean + Rocq proof … in ~3 minutes". I could not verify this claim at its source.

**Social media** [TECH/CO]
- @lanyon_ai's first post (2026-07-17 13:04 UTC) had 298 likes — [X](https://x.com/lanyon_ai/status/2078103512081117257)
- Gorard's quote-post (2026-07-17 13:14 UTC), "Now I can finally reveal why I've been quiet for so long! We're building a radically new kind of formally verified AI for science, math, engineering, and everything else, over at @lanyon_ai", had 1,136 likes and 36 replies — [X](https://x.com/getjonwithit/status/2078105967716082089)
- Follower counts could not be retrieved; X syndication was rate-limited.

**Conferences**
- No Lanyon talks or booths were found as of 2026-10-04.
- The Dorland Prize will be "announced at the 2026 APS Division of Plasma Physics Meeting in Chicago". Lanyon is a donor, not a listed presenter [INST] — [APS](https://www.aps.org/about/support/support-dorland-prize-computational-plasma)

### Inferences
- **Headcount** is most likely about 3–5. There is no public evidence of any employee beyond the founders.
- **Revenue.** The company is almost certainly pre-revenue. Signs: Dimension's "those things will come", no named customers, and the founders' own "tall order" remark.
- **Hiring as a technology signal.** The roles name formal verification, computational physics, applied mathematics and AI infrastructure. This fits a DSL-plus-proof-generation stack that targets physics solvers.
- **Engineering concentration.** Commits come mostly from the CEO's account, which suggests the public engineering output is concentrated in one founder. That may reflect who publishes rather than who builds.

### Gaps
- True headcount and payroll; any LinkedIn company page under another slug.
- X and LinkedIn follower counts; podcast appearances specifically about Lanyon (none found).
- Any customer, pilot, LOI or revenue figure.

## 5. Timeline (founding, launches, pivots, announcements through Oct 2026)

### Takeaway
The path ran as follows:
- 2024 to Feb 2026: academic precursor work at Princeton and PPPL.
- 2026-04-04: GitHub org created.
- 2026-06-05: domain registered.
- 2026-07-17: stealth exit, with solver repos published the same day.
- July 2026: round closed.
- 2026-08-14: NJ foreign registration.
- 2026-08-17: funding press release.
- September 2026: blog posts on the Dorland Prize and the company's vision, plus repo updates.

No pivots and no announcements between 2026-10-01 and 2026-10-04 were found.

### Cited Findings
| Date | Event | Evidence | Source |
|---|---|---|---|
| 2017 (summer) | Gorard's Wolfram Summer Research Institute project, "Automated Theorem Proving for Equational Logic" | [INST] | [Wolfram](https://education.wolfram.com/summer-research-institute/alumni/2017/gorard/) |
| 2014–2020 | Juno's PhD at UMD (Dorland/TenBarge) on continuum Vlasov–Maxwell methods, in collaboration with Hakim | [self-CV][CO] | [CV](https://spacalum.rice.edu/bios/james_juno_cv.pdf); [Lanyon blog](https://lanyon.ai/blog/dorland/) |
| Jan 2022 | Juno becomes PPPL Staff Research Physicist | [self-CV] | [CV](https://spacalum.rice.edu/bios/james_juno_cv.pdf) |
| Feb 2024 | Gorard joins Princeton RSE group / Gkeyll team at PPPL (led by Hakim) | [INST] | [Princeton RC (archived)](https://web.archive.org/web/20260206171050/https://researchcomputing.princeton.edu/about/people-directory/jonathan-gorard) |
| 2024-10-03 | First Gorard–Hakim–Juno joint paper (tetrad-first numerical GR) | [TECH] | [arXiv:2410.02549](https://arxiv.org/abs/2410.02549) |
| 2025-03-18/19 | "Shock with Confidence" (Gorard & Hakim): Racket pipeline generating verified C solvers with formal proofs; X thread | [TECH][CO] | [arXiv:2503.13877](https://arxiv.org/abs/2503.13877); [X](https://x.com/getjonwithit/status/1902158541839856071) |
| 2026-02-16 | BEACONS paper (Gorard, Hakim, Juno): formally verified neural PDE solvers; the founders' last arXiv paper to date | [TECH] | [arXiv:2602.14853](https://arxiv.org/abs/2602.14853) |
| 2026-04-04 | GitHub org `lanyonai` ("Lanyon AI Inc") created | [TECH] | [GitHub API](https://api.github.com/orgs/lanyonai) |
| H1 2026 (date unknown) | Company incorporated; Siqi Chen is first investor "after we incorporated" | [CO] | [Lanyon blog](https://lanyon.ai/blog/fundraising/) |
| 2026-06-05 | lanyon.ai domain registered (InternetX, privacy proxy) | [TECH] | [RDAP](https://rdap.identitydigital.services/rdap/domain/lanyon.ai) |
| ~2026-06-11 | @lanyon_ai X account created (decoded from user ID) | [TECH] | [X](https://x.com/lanyon_ai) |
| 2026-07-14/15 | Website build: sitemap `lastmod` dates (sitemap still lists `localhost:1313` URLs) | [TECH] | [sitemap.xml](https://lanyon.ai/sitemap.xml) |
| 2026-07-17 | **Stealth exit**: @lanyon_ai first post 13:04 UTC; Gorard's post 13:14 UTC; blog post "Announcing Lanyon AI" (dated "July 2026"); AdvectionDiffusion repo created 03:53 UTC | [CO][TECH] | [X](https://x.com/lanyon_ai/status/2078103512081117257); [X](https://x.com/getjonwithit/status/2078105967716082089); [Blog](https://lanyon.ai/blog/welcome/) |
| 2026-07-20 to 07-29 | Maxwell, Vlasov, Burgers, GR-Maxwell and compressible Euler solver repos published | [TECH] | [GitHub](https://github.com/lanyonai) |
| July 2026 | Initial round closes (Dimension lead; Industrious; Siqi Chen) | [CO] | [Lanyon blog](https://lanyon.ai/blog/fundraising/) |
| 2026-07-21 | Dimension announces $800M Fund III (context) | [PRESS] | [TechCrunch](https://techcrunch.com/2026/07/21/dimension-capitals-800m-third-fund-shows-the-intersection-of-science-and-compute-is-booming/) |
| 2026-07-24 | First Wayback capture of lanyon.ai, including /company/ with all three founder headshots | [TECH] | [Wayback](https://web.archive.org/web/20260724181740/https://lanyon.ai/) |
| 2026-08-14 | "LANYON AI INC." registered in NJ as a foreign for-profit corporation (Pennington) | [REG] | [NJ search](https://www.njportal.com/DOR/BusinessNameSearch/Search/BusinessName) |
| 2026-08-17 03:00 ET | PR Newswire: "Lanyon AI Emerges from Stealth…", $10.6M | [CO] | [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html) |
| Aug 2026 | Blog "Introducing our Partners" | [CO] | [Lanyon blog](https://lanyon.ai/blog/fundraising/) |
| 2026-08-24 / 08-27 | Fork of Gkeyll into the org; IdealMHD repo published | [TECH] | [gkeyll fork](https://github.com/lanyonai/gkeyll); [IdealMHD](https://github.com/lanyonai/IdealMHD) |
| 2026-08-31 | Wayback captures new research notes (ideal MHD, GR-Maxwell, Euler, nonlinear benchmarking) | [TECH] | [Wayback CDX](https://web.archive.org/cdx/search/cdx?url=lanyon.ai&matchType=domain) |
| Sep 2026 | Blog "The William D. Dorland Prize" (Juno); Lanyon donates to the APS prize endowment | [CO][INST] | [Lanyon blog](https://lanyon.ai/blog/dorland/); [APS](https://www.aps.org/about/support/support-dorland-prize-computational-plasma) |
| 2026-09-09 | AdvectionDiffusion and BurgersEquation repos updated ("Added new formally verified C implementations, Lean proofs of correctness…") | [TECH] | [Commits](https://github.com/lanyonai/AdvectionDiffusion/commits) |
| By 2026-09-13 | Blog "Our Vision" (Gorard), first captured 2026-09-13 | [CO][TECH] | [Lanyon blog](https://lanyon.ai/blog/vision/) |
| 2026-09-23 | Last push to the gkeyll fork | [TECH] | [GitHub](https://github.com/lanyonai/gkeyll) |
| 2026-10-01 to 10-04 | No new announcements found | — | searches, 2026-10-04 |

### Inferences
- The "stealth" phase was short in public terms. GitHub org (April) → domain (June) → launch (July) is about 3.5 months, though the underlying research ran from 2024 or earlier.
- **Scope widening, not a pivot.** The July 2026 blog names aerospace, propulsion and nuclear. The August press release also lists "GPU kernel optimization, frontier AI inference", and Dimension's quote mentions "chip tape-outs". The "Our Vision" post (September) speaks of building "a type system for the physical universe".

### Gaps
- The exact incorporation date.
- The exact publication dates of the blog posts: the pages show month only.
- Whether benchmark results were "released in the coming weeks" as the July post promised. Two benchmarking notes exist on the site, but they are undated.

## 6. Red flags and open questions

### Takeaway
Nothing suggests fabrication. The founders are real, well-published researchers. The company's existence and funding are corroborated by an NJ registry entry, an independent APS donor listing, the investors' own profiles, and a dense public GitHub record.

The concerns are as follows:
- absolute marketing claims that the company itself narrows in a footnote
- self-published, unreplicated benchmarks against named frontier models
- an unsupported "award-winning" descriptor
- missing legal details (home state, incorporation date, no Form D)
- stale institutional pages that leave departure dates and IP status unclear
- target markets widened between the July blog and the August release
- a very young domain and accounts, with no customers named

### Cited Findings
**Absolute correctness claims that the company itself narrows**
- The welcome post (July 2026) says "it is mathematically impossible for Lanyon ever to make a mistake" and that Lanyon is "100% reliable" [CO] — [Lanyon blog](https://lanyon.ai/blog/welcome/)
- Its own footnote narrows this: "It is perhaps more accurate to say that it is mathematically impossible for Lanyon ever to commit a misformalization … Lanyon guarantees perfect syntactic correctness. The issue of semantic correctness (i.e. how to guarantee that the formal specification actually matches the user's natural language intent) remains an open area of research for us."
- The press release repeats a version of the claim ("it is mathematically impossible for Lanyon's generated code not to match its formal specification"; "100% reliable") without the caveat [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)

**Unreplicated comparative claims**
- The company claims "only a tiny fraction of the token and compute cost of frontier models like GPT-5.6 and Fable 5" [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)
- Lanyon wrote "We will be releasing public benchmark scores attesting to this fact in the coming weeks" [CO], July 2026 — [Lanyon blog](https://lanyon.ai/blog/welcome/)
- The only benchmarks found are self-published, undated research notes. No independent replication was found as of 2026-10-04 — [Research notes](https://lanyon.ai/research/)

**Unsupported or inflated descriptors**
- "Award-winning applied mathematician" (Gorard) [CO]. Searches on 2026-10-04 found no named award. His Princeton bio lists an MPhil and no doctorate — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html); [Princeton RC (archived)](https://web.archive.org/web/20260206171050/https://researchcomputing.princeton.edu/about/people-directory/jonathan-gorard)
- Puffery such as "world-leading Princeton mathematicians and physicists" and "decades of formal methods expertise" (Dimension quote) [CO], 2026-08-17 — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)

**Legal and filing opacity**
- NJ shows a foreign corporation registered 2026-08-14 [REG], but no home state or incorporation date — [NJ search](https://www.njportal.com/DOR/BusinessNameSearch/Search/BusinessName)
- No Form D on EDGAR as of 2026-10-04 [REG] — [EDGAR FTS](https://efts.sec.gov/LATEST/search-index?q=%22Lanyon%20AI%22&dateRange=custom&startdt=2026-01-01&enddt=2026-10-04)
- The only valuation figure ($42M post) comes from Dealroom [DB] — [Dealroom](https://dealroom.co/companies/lanyon-ai/)

**A very young footprint**
- Domain registered 2026-06-05 behind a privacy proxy [TECH] — [RDAP](https://rdap.identitydigital.services/rdap/domain/lanyon.ai)
- X account created about 2026-06-11 [TECH] — [@lanyon_ai](https://x.com/lanyon_ai)
- GitHub org created 2026-04-04 [TECH] — [GitHub API](https://api.github.com/orgs/lanyonai)
- The address is a Regus suite — [Regus](https://www.regus.com/en/us/433)
- Minor site-hygiene signs [TECH] — [sitemap.xml](https://lanyon.ai/sitemap.xml); [Wayback CDX](https://web.archive.org/cdx/search/cdx?url=lanyon.ai&matchType=domain):
  - The sitemap lists `http://localhost:1313/` URLs.
  - Wayback captured Hugo dev-server `livereload.js` requests.
  - A linked `/code/ideal-mhd.mac` file returns 404.

**Unclear separation from Princeton/PPPL**
- On 2026-08-02, PPPL still listed Hakim as "Principal Research Physicist" and Deputy Head [INST] — [PPPL (archived)](https://web.archive.org/web/20260802195838/https://www.pppl.gov/people/ammar-hakim)
- On 2026-10-04, Hakim's GitHub bio still says "Lead Physicist for the Gkeyll Group" [TECH] — [GitHub](https://github.com/ammarhakim)
- On 2026-10-04, Juno's GitHub and Scholar still show PPPL [TECH] — [GitHub](https://github.com/JunoRavin)
- The precursor verification work, "Shock with Confidence" (2025) and BEACONS (2026), was done while the founders were Princeton/PPPL staff [TECH] — [arXiv:2503.13877](https://arxiv.org/abs/2503.13877); [arXiv:2602.14853](https://arxiv.org/abs/2602.14853)
- Lanyon forked the open-source Gkeyll repository (MIT-licensed) on 2026-08-24 [TECH] — [gkeyll fork](https://github.com/lanyonai/gkeyll)

**Target markets widened**
- The July blog's initial targets are "aerospace engineering, space and atmospheric propulsion, nuclear energy" [CO] — [Lanyon blog](https://lanyon.ai/blog/welcome/)
- The August press release adds "GPU kernel optimization, frontier AI inference" [CO] — [PR Newswire](https://www.prnewswire.com/news-releases/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing-302852383.html)

**Open outputs, closed generator (not a red flag; for context)**
- The generated solver code and proofs are public, and mostly MIT-licensed; only AdvectionDiffusion has no licence file [TECH], 2026-10-04 — [github.com/lanyonai](https://github.com/lanyonai)
- The DSL, compiler and agent that produce them are described as proprietary and are not published [CO] — [Lanyon blog](https://lanyon.ai/blog/welcome/)
- So outsiders can inspect and re-check the generated proofs, but cannot reproduce the generation step.

**No independent press or customers**
- Coverage consists of the company's own release and its rewrites [PRESS] — [TheSaaSNews](https://www.thesaasnews.com/news/lanyon-ai-emerges-from-stealth-to-build-the-future-of-scientific-and-technical-computing/); [TechEdge AI](https://techedgeai.com/lanyon-ai-emerges-with-10-6m-bet-on-provably-correct-scientific-ai/)
- No customers are named (see §4).

### Inferences
- **Overall credibility.** I see no evidence that the company is a façade. The research lineage (2025–26 arXiv papers) closely matches the product, and the GitHub artifacts are substantive.
  - The main credibility risk is the gap between the absolute marketing language and the narrower guarantee (syntactic correctness: implementation matches specification).
  - The question of whether a specification is the right one — the semantic gap — is openly unsolved by the company's own admission.
- **IP and separation.** Questions here are open, not red flags. It is unknown whether Princeton or PPPL holds rights in the pre-company pipeline, such as the Racket verification pipeline from "Shock with Confidence", and whether any founder kept a part-time or affiliate appointment. Neither is disclosed.
- **The missing Form D** is worth watching. Its continued absence after a disclosed $10.6M close would be unusual but has benign explanations; see §3.

### Gaps
- Any third-party technical validation, customer testimonial or replication of the benchmarks.
- The state and date of incorporation; any Form D; board composition.
- Founders' formal separation dates from Princeton and PPPL, and any licence or IP agreement.
- What "award-winning" refers to.
- X and LinkedIn audience metrics; whether any conference talk is scheduled (none found).
