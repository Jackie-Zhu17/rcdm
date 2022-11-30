# rcdm
Random Covariance Dynamic Model

## Description

This R package implements a random covariance dynamic model (RCDM) for joint estimation of sparse time-varying precision matrices.

</p>

<h2 id="overview-main">

Overview of main functions

</h2>

<table>
<colgroup>
<col style="width: 28%" />
<col style="width: 71%" />
</colgroup>
<thead>
<tr class="header">
<th>Function</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="even">
<td><code>randDynCov</code></td>
<td>Implements the Random Covariance Dynamic Model (RCDM) for joint estimation of sparse time-varying precision matrices. Optimization is conducted using block coordinate descent.</td>
</tr>
<tr class="odd">
<td><code>randDynCov_cluster</code></td>
<td>Implements the Random Covariance Dynamic Model (RCDM) for joint estimation of sparse time-varying precision matrices with pre-defined clusters(knots). </td>
</tr>
<tr class="even">
<td><code>bic_cal_m</code></td>
<td>Calculates the modified BIC for the RCDM.</td>
</tr>
</tbody>
</table>

<h2 id="install">

Installation

</h2>

Install the latest version of the package from GitHub with the following R code:

    if("devtools" %in% installed.packages() == FALSE) {
        install.packages("devtools")
    }
    
    devtools::install_github("Jackie-Zhu17/rcdm")

<h2 id="examples">

Example

</h2>

For a detailed simulation example for implementing the Random Covariance Dynamic Model (**RCDM**), see this.
