# halfrule (Python)

Optimal eligibility criteria for randomised trials. See the [main README](https://github.com/Marwan-Alqumbaey/halfrule) for the rules and examples.

```bash
pip install "git+https://github.com/Marwan-Alqumbaey/halfrule#subdirectory=python"
```

```python
from halfrule import eligibility

fit = eligibility(effect=[1.0, 0.6, 0.4, 0.1], weight=[0.2, 0.3, 0.3, 0.2])
fit.included
# array([ True,  True,  True, False])
```

Run the tests with `pytest`. The file `tests/fixtures.json` holds results of the R package, written by `tests/fixtures.R` (which needs the R package `jsonlite`), and the tests check that both packages agree.
