# Challenge_2
**Repository for Scripts and simulations for the Dream Challenge**

This repository contains the scripts used to generate the simulations for the SubChallenge 2 of the Allen Institute Cell Lineage Reconstruction Dream Challenge.

## About the main script

The matlab script "Sims_complete_lineage_Intertarget.m" is the one that generates the test dataset used for the Challenge.
It simulates the accumulation of mutations in a CRISPR cassete of 200 targets on a L3 larval stage hermaphrodite *C. elegans* (SubChallenge 2).

The script needs to be executed with matlab from its folder. For running from the terminal (in Linux) the following command should be used (matlab version R2019a or newer):

`matlab -batch "Sims_complete_lineage_Intertarget"`

The script will use the cell lineage tree file of C elegans with 1092 cells located in ./training_set/1000_cells.json

The main output will be the file FULLTREE_200_targets_30_states_0001rep_Dropout.fas in the folder "./simulations/" (folder needs to be created before), which will contain the cell ID and barcode of all the cells, like this:

> c_0000	0D0A0AA0G0A00A0000AA00A0A00A00E0AB00BC00AA00A0B0AEAAAA00000000A000000A000AI------A00000000A000AA0000A----------AD0000AAA00A0000A000AA0A0A00000EB00A00A0DA000000AA0A00AC0000A00000000D000A000A00C0A000000<br />
c_0001	A00B0A00ACAB0D0000A---------------A00A00GA0A0000000AA000000C000A0000AA00AA00A0AA0BAA000000000AE000000000A0AA0A00000000AA00AA0A0AA00AA00000C00A0B00A00A000000BA00A00B0AACACA0C-----------A00000B00EA000A0<br />
....<br />
c_1091	AEAA000AAFA0C--AABAAH0C0A00ALAAGAA0AAF0A0BE0DAAA00AFAG00CA0A0EC0B0A0AGB0A0A00A0A0AA00ABABAEABAABEAAB0ADAA---C0FA0EAAAADABAABFAEAAAACB0C0HCABAABAAAABAA00AAB0AAA00A0A0CFBAA00AAADAA0AAAA00AAAAAAAAA00AAAB<br />

for the challenge we took the first 1000 cells (until c_0999).

## Coding mutated/unmutated states

The mutated/unmutated states of each site were coded in a character state matrix, with each of the possible 30 mutations represented by letters "A-Z" and "a-d", following a random gamma distribution. If the target was missing due to an inter-target deletion event or a dropout, the character becomes "-". Once a target is mutated, it can no longer change, either to revert to the unmutated state or to transit to a new state.
 
The mutation events were implemented as following a Poisson distribution, with a probability of ~3x10-4 mutations per minute. Under the Poisson model, given a mutation rate μt (per unit time), the probability that a site remains unmutated after t minutes is:  e-(𝜇t.t). The general implementation of these simulations is similar to the ones described in (Salvador-Martínez et al., 2019).

## Inter-target deletions
 
In *C. elegans*, if two mutations occur in close targets, (no more than 20 targets apart in the recording array), within a short interval of time (~3 min) during a given cell division, all the targets between them are removed.

## Variables 

To modify the main parameters in the script the relevant variables can be modified in the lines 7-19 in the script (as shown in here)

> % Define the MAIN parameters<br />
%-------------------------------------------------------------------------<br />
% Number of targets<br />
targets = 200;<br />
% Mutation rate per minute<br />
Lambda = 0.00015;<br />
Norm_time = 1.88; % this is a normalising factor to get the real mins
%NUmber of sims<br />
Nsims = 100;<br />
% the number of states in the analysis. PAUP alows max 64 states.<br />
states = 30;<br />
% max_time<br />
max_time = 3065;<br />


**NOTE**
The cell lineage tree of the L3 larval stage cell lineage of the hermaphrodite C. elegans, used for this simulation comes from http://wormweb.org/celllineage as a json format (also available in this repository). In the annotated tree the cell division time was originally scaled by a factor of 1.88, so the real time of a given cell division can be obtained by dividing the time reported in the json file by 1.88. This is done with the variable "Norm_time". 
The simulation starts with one cell at minute 0 (fertilization) and ends with 1092 cells at minute 3,065 (L3 stage at 20ºC), as cell death was not considered.

This explains why the Mutation rate per minute in this script is 1.5x10-4, as this divided by 1.88 gives ~3x10-4, as reported above. Any desired modification to the Mutation rate will need to consider this.
