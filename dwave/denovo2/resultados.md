En este segundo experimento en D-Wave hemos hecho lo siguiente:

- **Crear la reconstrucción del gen a partir del resultado del DWave**: Es decir, coser de vuelta las lecturas en el orden dado por DWave para poder comparar la solución final con la inicial.

- **Automatizar la creación de tests**: Para ello creamos una cadena aleatoria con las letras 'ACGT' y la cortamos con tamaño de lectura fija y overlap entre lecturas fijo. Obtenemos así el gen inicial y las lecturas. Esto ha sido debidamente testeado.

- **Automatizar la realización de tests**: A partir del test creado anteriormente se realizan 10 tests de forma con tamaño total de gen 19, tamaño de cada lectura 10 y overlap entre lecturas 7. Estos valores son los mismos que en el ejemplo del paper QuASer. Esto resulta en 4 lecturas y un nodo de 4^2=16 elementos.

Los resultados obtenidos han sido los siguiente:

```
Test 1 results: Sucess
Test 2 results: Sucess
Test 3 results: Sucess
Test 4 results: Sucess
Test 5 results: Sucess
Test 6 results: Sucess
Test 7 results: Sucess
Test 8 results: Sucess
Test 9 results: Sucess
Test 10 results: Sucess
```

Además, hemos medido el tiempo de ejecución medio de cada test. Esto incluye unicamente la resolución del problema: A partir de las lecturas, obtener el gen final recreado (excluyendo la creación del test). El tiempo medio por test ha sido de `4.279282665252685` segundos, ejecutado en mi máquina personal.

Éxitos del experimento:
- Poder simular DWave localmente sin variar siquiera la API.
- Automatizar la creación y ejecución de tests.