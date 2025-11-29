# **stah — Survival analysis using average hazard with survival weight in Stata**

[![DOI](https://zenodo.org/badge/1106078435.svg)](https://doi.org/10.5281/zenodo.17756397)

**Authors:** Emily Xing, Lu Tian, and Hajime Uno

**License:** MIT License

---

## **Overview**

`stah` is a Stata command for performing **Average Hazard with Survival Weight (AH)** analysis, a robust and interpretable alternative to hazard‐ratio–based survival comparisons.
It implements AH-based methods for (a) one-sample estimation and (b) two-sample comparisons (unstratified and stratified). 

This repository includes:

* `ado/stah.ado` – the main Stata command
* `ado/stah.sthlp` – the help file
* `stah.pkg` – package definition for installation
* `stata.toc` – table of contents for the package site
* `example/` – example do-file

---

## **Installation**

You can install the package directly in Stata using:

```stata
net from "https://uno1lab.com/stata-stah/"
net install stah
```

After installation, view the documentation using:

```stata
help stah
```

---

## **Quick Example**

```stata
* Load PBC data
clear all
import delimited "pbc_data.csv", clear

* rescale time from day to year
gen year = time/365.25

* Create event variable (use "2" only)
gen event = (status == 2)

* Convert string treatment to numeric
gen trt_num = .
replace trt_num = 1 if trt == "1"
replace trt_num = 2 if trt == "2"

// Set survival data
stset year, failure(event)

// Two-sample analysis (treatment comparison)
stah trt_num, tau(5) reference(2)
```

This runs AH comparison up to τ = 5 years using the trt_num variable. More examples are provided in `example/example-stah.do`.

---

## **License**

This software is released under the **MIT License**.

© 2025

* Emily Xing
* Lu Tian
* Hajime Uno
