# PTM-Placenta-Atlas

! replace w/ paper link
! replace w/ GEO access

## Overview
This repository contains the code and instructions to reproduce the results from:
> **"A single-cell transcriptomic atlas of the pigtail macaque placenta in late gestation "**, Amanda Li 1,2, Richard Li3,, Hazel Huang1,, Hong Zhao1, Briana Del Rosario1, Miranda Li1,4, Edmunda Li1, Andrew Vo1,4, Gygeria Manuel1,5, Orlando Cervantes1,6, Raj Kapur7,8, Jeff Munson9, Austyn Scanlon10, Michelle Coleman10, Melissa Berg11, Britni Coman11, Brenna Menz11, Jin Dai11, Inah Golez11, Solomon Wangari11, Chris English11, Audrey Baldessari11, Lakshmi Rajagopal6,10, John Cornelius1,\*, Kristina Adams Waldorf 1,6,10,\* 

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