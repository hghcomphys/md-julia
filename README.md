# Molecular Dynamics in Julia

In this repository, I've crafted a simple molecular dynamics code in Julia that simulate a Lennard-Jones interacting atoms at room temperature. This project represented a valuable chance for me to delve deeper into this high-level but fast programming language, and by sharing I hope it proves to be beneficial for you as well.

## Example

The `md.jl` file within the `examples` directory demonstrates the simulation of 100 Helium atoms within the simulation box.

To execute this example, please utilize the following command:
```bash
$ julia md.jl

Molecular Dynamics in Julia
===========================
Potential: Lennard-Jones
Precision: Float64
Number of atoms: 100
Temperature: 300.0
Integration: NVT
0          t[ps]:0.00000    Temp[K]:298.0642697858  Etot[Ha]:0.1414726142    Epot[Ha]:-0.0001143922   Pres[kb]:0.1905196184   
100        t[ps]:0.05000    Temp[K]:162.5451267265  Etot[Ha]:0.0771101440    Epot[Ha]:-0.0001023240   Pres[kb]:0.1039384608   
200        t[ps]:0.10000    Temp[K]:351.5004051327  Etot[Ha]:0.1668776566    Epot[Ha]:-0.0000926755   Pres[kb]:0.2246785962   
300        t[ps]:0.15000    Temp[K]:398.1503094939  Etot[Ha]:0.1890531454    Epot[Ha]:-0.0000769052   Pres[kb]:0.2544489953   
400        t[ps]:0.20000    Temp[K]:306.6724819343  Etot[Ha]:0.1456073578    Epot[Ha]:-0.0000687365   Pres[kb]:0.1959680985   
500        t[ps]:0.25000    Temp[K]:315.3815996065  Etot[Ha]:0.1497436227    Epot[Ha]:-0.0000694917   Pres[kb]:0.2015429553   
600        t[ps]:0.30000    Temp[K]:296.7910968391  Etot[Ha]:0.1409128136    Epot[Ha]:-0.0000694080   Pres[kb]:0.1896902351   
700        t[ps]:0.35000    Temp[K]:305.4121319993  Etot[Ha]:0.1450084830    Epot[Ha]:-0.0000689176   Pres[kb]:0.1952454235   
800        t[ps]:0.40000    Temp[K]:273.6631686965  Etot[Ha]:0.1299249339    Epot[Ha]:-0.0000710191   Pres[kb]:0.1749121103   
900        t[ps]:0.45000    Temp[K]:294.5114240205  Etot[Ha]:0.1398254889    Epot[Ha]:-0.0000738386   Pres[kb]:0.1881603751   
1000       t[ps]:0.50000    Temp[K]:318.0675487984  Etot[Ha]:0.1510167955    Epot[Ha]:-0.0000722032   Pres[kb]:0.2032178801   
```