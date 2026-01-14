# Fix: GOD-Berechnung korrigieren

## Problem
Die Kurven im Baur-Glässner-Diagramm verlaufen in die falsche Richtung.

## Ursache
Die GOD-Formel ist falsch. Aktuell:
```python
GOD = K / (1 + K)  # FALSCH
```

## Lösung
Ändere in `thermodynamics.py` die Funktion `GOD_from_K`:

```python
def GOD_from_K(K, stoichiometry=1):
    """
    K = (p_oxidized / p_reduced)^n
    GOD = 1 / (1 + K^(1/n))
    """
    K_eff = K ** (1.0 / stoichiometry)
    return 1.0 / (1.0 + K_eff)  # <-- Invertiert!
```

Auch in der Methode `GOD()` der Klasse `BaurGlassnerThermo` anpassen:

```python
def GOD(self, reaction_name: str, T_kelvin: float) -> float:
    # ... (Gültigkeitsprüfung bleibt gleich)
    
    K = self.equilibrium_K(reaction_name, T_kelvin)
    stoich = rxn['stoich']
    K_eff = K ** (1.0 / stoich)
    return np.clip(1.0 / (1.0 + K_eff), 0, 1)  # <-- Hier ändern!
```

## Erklärung
- K > 1 bedeutet: Reduktion ist begünstigt → niedriger GOD (reduzierende Atmosphäre)
- K < 1 bedeutet: Oxidation ist begünstigt → hoher GOD (oxidierende Atmosphäre)
