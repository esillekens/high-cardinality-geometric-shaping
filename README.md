# High-Cardinality Geometrical Constellation Shaping for the Nonlinear Fibre Channel
[![Generic badge](https://img.shields.io/badge/doi:-10.1109/JLT.2022.3197366-blue.svg)](https://doi.org/10.1109/JLT.2022.3197366) 
### Eric Sillekens , Gabriele Liga , Domanic Lavery , Polina Bayvel , Robert I. Killey (2022)
----

Code used to generate the results for the paper titled "High-Cardinality Geometrical Constellation Shaping for the Nonlinear Fibre Channel". Contains the GMI calculation including the gradient, the functions for the channel-specific chain rule and an implementation of the trust-region algorithm with SR1 Hessian update.



## Usage

Optimise a 2D constellation for the AWGN channel
```matlab
optimiseAWGN.m
```

Optimise a 4D constellation for the AWGN channel
```matlab
optimise4DAWGN.m
```

Optimise a 2D constellation for the nonlinear fibre channel
```matlab
optimiseAWGN.m
```