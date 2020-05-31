# SSEIH model

'SSEIH.R' contains the main code for running the model. It sources 'betasSSEIH.R' and 'inputsSSEIH.R'. 
'plotSSEIH.R' was used to create the plots for the [report published by SSI on May 6th 2020](https://files.ssi.dk/Ekspertrapport-af-den-6-maj).

### File "data_for_kaare2.csv" not available - sample file is now provided
The data for hospitalizations by region and age that is loaded as 'mdt' cannot be shared at this point. Sharing is pending legal clarification.
Meanwhile a file with the same structure is shared as "data_for_kaare_na.csv". All numbers are replaced by NAs.

### Updated to cover report on May 20th 2020
The code is updated to the version that was used for the report on May 20th 2020.
The main update is that extra phases and summer break are added.
The parameter 'p' was treated as a vector and it is converted to a list 'pl' to actually use the values in 'pl' rather than in the surrounding scope.
The code for making the beta matrices is provided in the Contacts folder