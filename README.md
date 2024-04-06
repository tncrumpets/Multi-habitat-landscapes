This repository contains R-scripts to run the analyses associated with the following manuscript:

"Multi-habitat landscapes support diversity, stability and improved function" by Hackett T.D., Sauve A.M.C, Maia K.P., Montoya D., Tylianakis J., Potts S., Vaughan I., Memmott J.


Analyses are separated in the following folders:

* `emergent_properties` contains R scripts to test whether the network properties we measured
on sites combining multiple habitats are different than expected if landscape scale food
webs were only the sum of their habitats.

* `robustness_analysis` contains R scripts to calculate the robustness of landscape-scale food webs to 1) the extinction of the least to most common plant species (landscape homogeneisation), and to 2) random extinctions (species loss *per se*).

* `sampling_completeness` gathers R scripts to calculate sampling completeness in the sampled networks, and their effects on the results presented in the manuscript.

All of these folders are organised the same way : a `main_scripts` folder listing one script per analysis and an `outputs` folder for any figures or table resulting for the analysis. In some cases, a `functions` folder gathers each function that can be called by the different R-scripts in the `main_scripts` folder.

The data analysed in this manuscript are gathered in the folder `data`.


