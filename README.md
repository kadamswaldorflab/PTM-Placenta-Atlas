# PTM-Placenta-Atlas

! https://www.biorxiv.org/content/10.1101/2025.08.15.669966v1
! replace w/ GEO access

## Overview
This repository contains the code and instructions to reproduce the results from:
> **"A single-cell transcriptomic atlas of the pigtail macaque placenta in late gestation "**, Amanda Li, Richard Li, Hazel Huang, Hong Zhao, Briana Del Rosario, Miranda Li, Edmunda Li, Andrew Vo, Gygeria Manuel, Orlando Cervantes1,6, Raj Kapur7,8, Jeff Munson9, Austyn Orvis, Michelle Coleman, Melissa Berg11, Britni Curtis, Brenna Menz, Jin Dai, Inah Golez, Solomon Wangari, Chris English, Audrey Baldessari, Lakshmi Rajagopal, John Cornelius, Kristina Adams Waldorf 

## Repository Structure
- `src/` — Source code.
- `results/` — Output figures/tables from the paper.
- `docs/` — Additional documentation.

## Installation (R packages)
```bash
git clone https://github.com/kadamswaldorflab/PTM-Placenta-Atlas.git
cd PTM-Placenta-Atlas
```
```R
pkgs <- readLines("requirements-R.txt")
install.packages(pkgs)
```

## Installation (Python packages)
```bash
cd PTM-Placenta-Atlas
pip install -r requirements-python.txt
```
## USAGE
```bash
