{smcl}
{* 26November2025}{...}
{cmd:help stah} 
{hline}

{title:Title}

{p2colset 5 16 18 2}{...}
{p2col :{hi:stah} {hline 2}}Comparing survival distributions using average hazard{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 15 2}
{cmd:stah}
[{it:groupvar}]
{ifin}
[{cmd:,} {opt tau(#)} {opt strata(varname)} {opt level(#)} {opt reference(#)}]


{title:Description}

{pstd} {cmd:stah} estimates the average hazard (AH) as a summary measure of survival distributions over pre-specified time window [0,τ].
The AH can be interpreted as expected "events per person-time" up to truncation time τ, providing
a model-free, censoring-robust alternative to conventional hazard ratios.

{pstd} {cmd:stah} supports single-arm and two-sample analyses, reporting the AH estimate between groups with standard errors 
and confidence intervals. For two-sample comparisons, specify {it:groupvar} (a numeric variable such as 0/1) to obtain group-specific 
AH estimates and between-group contrasts: the difference in AH (DAH) and the ratio of AH (RAH).

{pstd} {cmd:stah} also supports stratified analysis via direct standardization using the
{opt strata()} option, thus providing adjusted estimates without assuming homogeneous treatment
effects across strata. The command presumes the data have first been declared as survival-time
data using the {cmd:stset} command.


{title:Options}

{phang}{opt tau(#)} A scalar value to specify the truncation time point for the AH
calculation. Tau must be less than or equal to the largest observed time (either event or censor) in each treatment group. 
The default value is the minimum of the largest observed event time across all treatment groups.

{phang}{opt strata(varname)} Performs stratified analysis using the specified stratification
variable. The variable takes integer values representing different strata. When specified, {cmd:stah} conducts direct standardization
with proportional stratum weights proportional to observed stratum sizes; otherwise, the default is proportional to the stratum size
in the input data. 

{phang}{opt level(#)} set confidence level; the default is level(95).

{phang}{opt reference(#)} specifies the reference group in two-sample analysis. The default
is the group with the smallest {it:groupvar} value.

{phang} {opt by} is allowed with {cmd:stah}; see {manhelp by D}.


{title:Examples}

{p 8 14 2}{cmd:. stah, tau(5)}

{p 8 14 2}{cmd:. stah treatment, tau(5)} 

{p 8 14 2}{cmd:. stah treatment, tau(5) reference(2)} 

{p 8 14 2}{cmd:. stah treatment, tau(5)  reference(2) strata(bili_strata)}

{p 8 14 2}{cmd:. stah treatment, tau(10) strata(stage) level(90)}


{title:Authors}

{pstd}Emily Xing{p_end}
{pstd}Harvard University{p_end}
{pstd}exing at college.harvard.edu{p_end}

{pstd}Lu Tian Uno{p_end}
{pstd}Stanford University{p_end}
{pstd}lutian at stanford.edu{p_end}

{pstd}Hajime Uno{p_end}
{pstd}Dana-Farber Cancer Institute{p_end}
{pstd}huno at ds.dfci.harvard.edu{p_end}

{title:Also see}

{p 4 14 2}Uno H, Horiguchi M. Ratio and difference of average hazard with survival weight: 
New measures to quantify survival benefit of new therapy. Stat Med. 2023 Jan 5;42(7):936–52.

{p 4 14 2}Qian Z, Tian L, Horiguchi M, Uno H. A novel stratified analysis method for 
testing and estimating overall treatment effects on time-to-event outcomes using average 
hazard with survival weight. Statistics in Medicine 2025, 44(7), e70056.

{p 4 14 2}R package {cmd:survAH}: Uno H, Horiguchi M, Qian Z. uno1lab/survAH; 2025. Version 1.1.2. Available from: {browse "https://doi.org/10.5281/zenodo.15667795"}


