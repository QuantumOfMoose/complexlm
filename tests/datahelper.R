###
### Just a quick script to clean and package the Hall effect data as a .rda file
### William Ryan
### 19 June 2023

skip <- 19 # How many rows to skip before reading data from the csv file.
filename <- "data/12April2023_Copper_RT_vacrawdata.csv"
CuHallData <- read.csv(filename, header = T, row.names = NULL, sep = ",", skip = skip)  # Read in the Hall effect data.
dropcols <- grepl("Rot|V..Vrms|Vx.Noise|Theta|I..C|V..contact.1", colnames(CuHallData))
CuHallData <- CuHallData[,!dropcols] # Drop extraneous columns.
CuHallData$V..contact <- vapply(CuHallData$V..contact, FUN = function(x) {ifelse(x == 4, 'D', 'F')}, FUN.VALUE = 'A') # Replace the contact numbers with a letter indicating the contact arrangement
names(CuHallData)[4] <- "Contact.Arrangement" # And change the name of that column.
CuHallData$Vx..Vrms. <- sqrt(2) * complex(real = CuHallData$Vx..Vrms., imaginary = CuHallData$Vy..Vrms.) # Collect the voltage outputs into a single complex number, also convert from rms to amplitude.
names(CuHallData)[9] <- "OutputV"
CuHallData <- CuHallData[-10] # Ditch the now superfluous Vx..Vrms
CuHallData$RMS.Magnetic.Field..T. <- sqrt(2) * CuHallData$RMS.Magnetic.Field..T. # Convert magnetic field from rms to amplitude.
names(CuHallData)[3] <- "Magnetic.Field.T"
names(CuHallData) <- gsub("(.)\\1+", "\\1", names(CuHallData)) # Remove the repeated periods from the names.
CuHallData$Resistivity.Ohm.cm <- rep(2.54638e-06, length(CuHallData$Contact.Arangement)) # Manualy add the resistivity value calculated by the Van Der Pauw method from DC current-voltage sweeps.

### That aught to be sufficient cleaning. Now we can save this dataframe as a .rda file

save(CuHallData, file = "data/Copper_AC_Hall_Effect.rda")
