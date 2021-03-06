# FOLD : Fully powered method for OverLapping Data.
**FOLD/FOLD-split** is a free, open-source meta-analysis framework, which achieves optimal performance when case-control studies include shared control samples.

## Motivation
Recent studies have developed several meta-analysis methods that can systematically account for overlapping subjects at the summary statistics level.

Existing methods : Bhattacharjee et al. 2012; Bulik-Sullivan et al. 2015; Han et al. 2016; Lin and Sullivan 2009; Zaykin and Kozbur 2010

We observed that when not all controls are shared between studies, the existing approaches show a severe power drop compared with splitting(a strategy that splits genotype data of overlapping samples into individual studies).

![이미지](/image/power.jpg)

## Downloading the package at once
In order to download `FOLD`, you can clone this repository via the commands.

```
$ git clone https://github.com/hanlab-SNU/FOLD.git
$ cd FOLD
```


## Striking phenomenon

* Simulation code:

  [simulate-v1.0 (R script)](/FOLD_simulate/simulate-v1.0.R) ← Right click and download

* Command
  ```
  > source(“simulate-v1.0.R”)
  > simulate()
  ```

Inspired by this phenomenon, we propose a power-preserving meta-analysis framework for overlapping samples called **FOLD** (**F**ully-powered method for **O**ver**L**apping **D**ata).

## FOLD
FOLD obtains multiple summary statistics from each study, such that each statistic is calculated using samples that are homogeneous in terms of their information.

The frame of the gold standard approach and FOLD are

<img src="/image/FOLD_scheme.png" width="700">

### Code and manual Links:

[FOLD](/FOLD_data/FOLD.zip) ← Right click and download

## FOLD-split
FOLD can be complicated if the number of configurations of sharing(T) is large. In such situations, splitting can be simpler.

FOLD-split helps with split individuals to maximize power.

Equal splitting and Case-based splitting are two widely known methods, we compared FOLD-split with Equal and Case-based splitting approaches

<img src="/image/FOLD-split.jpg" width="700">

* Power of different methods:

Strategy A is a strategy that equally distributes shared controls into studies(Equal splitting)

Strategy B is a strategy that distributes shared controls proportionally to the number of cases in each study(Case-based splitting)

### Code and manual Links:

[FOLD-split](/FOLD_data/FOLD-split.zip) ← Right click and download

## Reference
[Kim EE, Lee S, Lee CH, Oh H, Song K, Han B. FOLD: a method to optimize power in meta-analysis of genetic association studies with overlapping subjects. Bioinformatics. 2017 Dec 15;33(24):3947-3954. doi: 10.1093/bioinformatics/btx463. PMID: 29036405; PMCID: PMC5860085.](https://academic.oup.com/bioinformatics/article/33/24/3947/3980249)

## Contact
Eunji Kim ([eunkim0116@gmail.com](mailto:eunkim0116@gmail.com))
