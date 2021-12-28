# Magi_TamR_2021

## Mathematical modeling of the TAM resistance acquisition process
Magi, S., Ki, S., Ukai, M. et al. A combination approach of pseudotime analysis and mathematical modeling for understanding drug-resistant mechanisms. *Sci. Rep.* **11**, 18511 (2021). https://doi.org/10.1038/s41598-021-97887-z

This model was constructed by using [BioMASS](https://github.com/biomass-dev/biomass) platform.

## Simulate model

```python
from biomass.models import Magi_TamR_2021 as model
from biomass import run_simulation

run_simulation(model, viz_type='average', show_all=False, stdev=True)
```
