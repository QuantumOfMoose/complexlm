#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom pracma geo_median
#' @import MASS
#' @import stats
#' @import mathjaxr
#' @importFrom graphics abline axis frame legend lines mtext panel.smooth par points strheight strwidth text title
#' @importFrom grDevices as.graphicsAnnot dev.flush dev.hold dev.interactive devAskNewPage extendrange
## usethis namespace: end
NULL

#' AC (Complex) Hall effect data measured from a thin-film copper sample.
#' 
#' \loadmathjax
#' This data serves as a 'real-world' example of complex valued measurements that can be analyzed by the methods contained 
#' in 'complexlm'. The Hall effect is a perpendicular (to the magnetic field and current) voltage that appears across an
#' electrically conductive material when it is placed in an external magnetic field. It is commonly used to determine the 
#' mobility and density of charge carriers (electrons and/or holes) within materials. 
#' 
#' The AC Hall effect is a more advanced technique that uses a time varying (sinusoidal) magnetic field and lock-in voltage 
#' measurement. The measured output signal is thus a wave, and best described by complex numbers. This strategy drastically 
#' lowers the noise floor, enabling measurement of difficult, low-mobility, materials. 
#' 
#' This data was take at room temperature in a vacuum chamber using a custom built and programmed system. A Keithley 2636
#' sourcemeter provided the excitation current, while an SRS 830 lock-in amplifier measured the voltage signal from the sample.
#' The sample consisted of a 0.8 mm square of 1000 Angstroms thick evaporated copper film on a 1.5 centimeter square die of oxidized crystalline silicon.
#' 
#' The variables included in the dataframe are as follows:
#' 
#' \itemize{
#'    \item{Temperature.K.}{The temperature of the sample during this measurement. Units of Kelvin.}
#'    \item{Magnet.Field.Frequency.Hz.}{The frequency of the oscillating magnetic field used for this measurement. Units are Hertz.}
#'    \item{Magnetic.Field.T.}{The amplitude of the magnetic field during this measurement. Units of Tesla.}
#'    \item{Contact.Arrangement}{This measurement involves four electrical contacts (current in and out, and voltage hi and lo), placed at the corners. 
#'          Corresponding contacts must be opposite each other, so there are two possible arrangements: "D" and "F".}
#'    \item{Current.Set.A.}{The desired current to be sent through the sample for this measurement. There was an error in recording beyond the 8th row. Units of Amperes.}
#'    \item{Current.In.meas.A.}{The current measured leaving the sourcemeter, towards the sample. Units of Amperes.}
#'    \item{Current.Out.meas.A.}{The current measured exiting the sample. Units of Amperes.}
#'    \item{Source.V.V.}{The voltage generated across the sample in order to produce the desired current. Units of Volts.}
#'    \item{OutputV}{The complex voltage measured across the sample. Ideally proportional to the current and magnetic field. Units of Volts.}
#'    \item{Current.leak}{The difference between the input and output currents. Units of Amperes.}
#'    \item{Resistivity.Ohm.cm}{The resistivity of the sample as calculated by the Van Der Pauw method from previous DC current-voltage sweeps. Units of Ohm*centimeter.}
#'    \item{thickness.cm.}{The thickness of the copper film in centimeters.}
#' }
#' 
#' @keywords datasets
#' @format A dataframe with 240 rows and 11 columns. Most names contain the units of the column after the last or second to last period.
#' @author William Ryan \email{wryan1@binghamton.edu}
#' 
"CuHallData"