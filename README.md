# PBN_pilot_analysis

End-to-end pilot analyses for the Psychiatric Biomarker Network (PBN): data intake, QC, and exploratory association workflows across WES and WGS cohorts.

> **Status:** active methods sandbox / pilot repo  
> **Data access:** controlled (PBN / IRB-approved only). No raw data are stored in this repository.

---

## Table of Contents
- [Overview](#overview)
- [Repo Structure](#repo-structure)
- [Getting Started](#getting-started)
- [Environment & Dependencies](#environment--dependencies)
- [Inputs & Expected Layout](#inputs--expected-layout)
- [Workflows](#workflows)
  - [WES](#wes)
  - [WGS](#wgs)
- [Results & Reproducibility](#results--reproducibility)
- [Project Notes](#project-notes)
- [Contributing](#contributing)
- [Citations](#citations)
- [License](#license)

---

## Overview

This repository hosts pilot pipelines and notebooks used to:
1. Validate intake/QC for PBN sequencing data (WES/WGS),
2. Prototype variant-level and sample-level QC thresholds,
3. Explore association signals and cohort harmonization,
4. Produce sharable figures/tables for internal reviews and PBN meetings.

Code is primarily in Jupyter notebooks (≈77%), with helper scripts in Shell (≈15%), R (≈7%), and light Python (≈1%). :contentReference[oaicite:1]{index=1}

---

## Repo Structure

