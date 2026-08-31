# About this repository
In this repository, I have included the code I used for my MGI thesis titled: Intergrating remote sensing and in-situ NO2 measurements to quantify dry nitrogen oxides deposition to veluwe forest.

The analysis pipeline is organized into seven stages. It begins with acquiring and pre-processing in-situ and TROPOMI NO₂ and meteorological data (A), followed by flux footprint modelling to define the spatial domain of representativeness (B) and its integration with remote sensing data (C). Vegetation indices derived from this integration are then analyzed (D, RQ1) alongside in-situ measurements (E, RQ1) to characterize surface conditions. These outputs feed into inferring surface NO₂/NOₓ concentrations (F, RQ2), which are finally used to estimate NO₂ dry deposition (G, RQ3).

In the Python and R, the project folder should be set as the working directory, and the environment should be activated first.

To download stacks from GEE, the resulted footprint area shapefile need to be imported as table, and defined as aoi.

There could be bugs and discrepancies due to the scale and limitation of time of the project, if you can not fix it, feel free to contact me!

Contact information: Zhiyu Wu, zhiyu.wu@wur.nl, zhiyuwu2023@gmail.com
