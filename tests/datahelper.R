###
### Just a quick script to clean and package the Hall effect data as a .rda file
### William Ryan
### 19 June 2023

skip <- 19 # How many rows to skip before reading data from the csv file.
filename <- "data/12April2023_Copper_RT_vacrawdata.csv"
rawdata <- read.csv(filename, header = T, row.names = NULL, sep = ",", skip = skip)  # Read in the Hall effect data.
dropcols <- grepl("Rot|V..Vrms|Vx.Noise|Theta|I..C|V..contact.1", colnames(rawdata))
rawdata <- rawdata[,!dropcols] # Drop extraneous columns.
rawdata$V..contact <- vapply(rawdata$V..contact, FUN = function(x) {ifelse(x == 4, 'D', 'F')}, FUN.VALUE = 'A') # Replace the contact numbers with a letter indicating the contact arrangement
names(rawdata)[4] <- "Contact.Arrangement" # And change the name of that column.
names(rawdata) <- gsub("(.)\\1+", "\\1", names(rawdata)) # Remove the repeated periods from the names.
rawdata$Resistivity.Ohm.cm <- rep(2.54638e-06, length(rawdata$Contact.Arangement)) # Manualy add the resistivity value calculated by the Van Der Pauw method from DC current-voltage sweeps.

### That aught to be sufficient cleaning. Now we can save this dataframe as a .rda file
CuHallData <- rawdata # First give it a better name.
save(CuHallData, file = "data/Copper_AC_Hall_Effect.rda")
