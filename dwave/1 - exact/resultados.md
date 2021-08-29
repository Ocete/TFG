En este primer experimento en D-Wave nos basamos en el ejemplo del paper QuASer leído, mirar allí los resultados esperados (página 12 tiene buenos diagramas).


Obtenemos los siguientes resultados:

```
{'n0t0': 0, 'n0t1': 0, 'n0t2': 1, 'n0t3': 0, 'n1t0': 0, 'n1t1': 0, 'n1t2': 0, 'n1t3': 1, 'n2t0': 1, 'n2t1': 0, 'n2t2': 0, 'n2t3': 0, 'n3t0': 0, 'n3t1': 1, 'n3t2': 0, 'n3t3': 0} --> -7.98113883008419
{'n0t0': 0, 'n0t1': 1, 'n0t2': 0, 'n0t3': 0, 'n1t0': 0, 'n1t1': 0, 'n1t2': 1, 'n1t3': 0, 'n2t0': 0, 'n2t1': 0, 'n2t2': 0, 'n2t3': 1, 'n3t0': 1, 'n3t1': 0, 'n3t2': 0, 'n3t3': 0} --> -7.98113883008419
{'n0t0': 1, 'n0t1': 0, 'n0t2': 0, 'n0t3': 0, 'n1t0': 0, 'n1t1': 1, 'n1t2': 0, 'n1t3': 0, 'n2t0': 0, 'n2t1': 0, 'n2t2': 1, 'n2t3': 0, 'n3t0': 0, 'n3t1': 0, 'n3t2': 0, 'n3t3': 1} --> -7.98113883008419
{'n0t0': 0, 'n0t1': 0, 'n0t2': 0, 'n0t3': 1, 'n1t0': 1, 'n1t1': 0, 'n1t2': 0, 'n1t3': 0, 'n2t0': 0, 'n2t1': 1, 'n2t2': 0, 'n2t3': 0, 'n3t0': 0, 'n3t1': 0, 'n3t2': 1, 'n3t3': 0} --> -7.98113883008419
```

**Recordemos que**: el índice `nitj` nos indica que la cadena (el nodo) `i` se coloca en la posición `j` (se recorre en el instante `j`) si solamente si vale 1.

Estas 4 son las soluciones encontradas por el sistema que tienen mínima energía, y todas tienen exactamente la misma energía. Estas corresponden a las secuencias:

- `2 -> 3 -> 0 -> 1`
- `4 -> 0 -> 1 -> 2`
- `0 -> 1 -> 2 -> 3`
- `1 -> 2 -> 3 -> 0`

Que son en realidad la misma solución, pues únicamente varía el nodo inicial. Éxito!
