# Generator for Heatsink Profiles
This project is a simple MVP on a method to dynamically generate profiles for surface-to-air heatsinks. It was supposed to be a few-day exercise in solving the heat-equation numerically and trying a simple approach to grow a heatsink. It was build in 2023.

It turned into a multi-week project due to a fatal choice early on: visualising the output in HTML-table to "speed up development".

## Theory of operation
The simulation implements the heat-equation at steady-state which reduces to computing the average for each cell and iterating until the heat-change between two iterations is below a chosen threshold.

Special care has to be taken at the aluminum-air-interface, where averaging has to be weighted by heat-conductance with each adjecent cell.

Since only the profile is generated, the problem can be solved in 2 dimensions. That approach is trivial for solid aluminum, but since air is moving through the heatsink, only laminar flow can be evaluated. Laminar flow can be solved using the heat-equation again, as each cell's velocity is the average of cells close by. Airflow is another source of heat-flow that has to be taken into account. This requires the length of the heatsink and design air-speed to be known. By knowing the air-speed in volume per second and the heat-flow in Joul per second, the change in heat-content of a single cell can be calculated.

After equillibrium is reached, the spot on the aluminum-surface with the highest temperature-difference to air gets another cell of aluminum attached.

## Usage
A rather untidy web-view in `heatsink.html` should load the `heatsink_simulator.js`. A colorful heatmap should be mostly black. By clicking "Grow on" the simulation is run to add several blobs of aluminum.

## Disclaimer
Although this project was done to implement physical principles, it was never tested whether it is physically accurate. There are known bugs (e.g. once the heatsink reaches any border) and no tests were written. If there is any interest in this project, I might rewrite this project into a more usable version.

## Word of advice
Write simple heat-equation solvers! It's simple, it's surprisingly versitile, and the difficult stuff is all edge-cases. (literally: the edges of the volume need very special care.)