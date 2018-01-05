# ellingham
*A Python script to generate Ellingham diagrams from energy of formation for binary compounds.*

By Levi Tegg

## Description
Ellingham diagrams plot the standard free energy of formation, \DeltaG<sub>f</sub><sup>o</sup>, as a function of temperature. For a description of the uses, I recommend *Free Energy of Formation of Binary Componds: An Atlas of Charts for High-Templerature Chemical Calculations* by Thomas B. Reed, published 1971 by MIT Press, Cambridge, Massachusetts.

The data in the chart has been obtained from the resources above, as well as *Thermodynamics of Binary Metallic Carbides: A Review* by R. G. Coltters in 1985, in *Materials Science Engineering* 76, pages 1-50.

## Requirements
The script was written for Python 3.5, but with minor adjustments it should work for Python 2.x. It uses the `numpy` and `matplotlib` libraries.

The .pdf and the .png can be downloaded and used as-is.

## Issues and Feedback
Any compound can be added to the diagrams if you know the start and end co-ordinate of the \DeltaG<sub>f</sub><sup>o</sup>(T) line.

In this instance, I have elected to embed the raw data in the script, rather than include it in a separate file. This is because the project started off small, then grew as I felt that more compounds should be added. In retrospect, this was unwise.
