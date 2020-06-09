# SSEIH model

'SSEIH.R' contains the main code for running the model. It sources 'betasSSEIH.R' and 'inputsSSEIH.R'. 
'plotSSEIH.R' was used to create the plots for the [report published by SSI on May 6th 2020](https://files.ssi.dk/Ekspertrapport-af-den-6-maj).

### Input file is now available - with data - after legal clarification. 
The file ’linelist_snapshot.csv’ is based on data from LPR3 and enriched with daily updated snapshots for COVID-19 patients in the regions. In the modelling the focus is on hospitalizations longer than 12 hours and it is therefore required that a patient is in snapshots from two consecutive days before that patient is counted as hospitalized with COVID-19. This excludes some patients with short visits (still longer than 12h) until they appear in LPR3.

### Updated to cover report on May 20th 2020
The code is updated to the version that was used for the report on May 20th 2020.
The main update is that extra phases and summer break are added.
The parameter 'p' was treated as a vector and it is converted to a list 'pl' to actually use the values in 'pl' rather than in the surrounding scope.
The code for making the beta matrices is provided in the Contacts folder