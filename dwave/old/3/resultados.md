# Test 1

**Ejecución completada en el D-Wave**: He conseguido lanzar experimentos en el D-Wave, aunque he descubierto que la creación de la solución es incorrecta, así que en el siguiente test corregiremos eso. Por ahora, cuál es le proceso de lanzamiento de un proceso en el D-Wave:

Preliminar: creamos el diccionario (solía ser `quboDict`, aqui es simplemente `Q`) que mapea lados del grafo con su energía como lo hacíamos antes.

Conectamos con la plataforma Leap para obtener el `solver` (arquitectura del grafo donde se lanzará el algoritmo) y el `sampler`, clase que utilizamos para resolver el problema. Para ello utilizamos un archivo de configuración muy simple que contiene nuestra api key:

```
config_file='../dwave.conf'
client = Client.from_config(config_file, profile='ocete')
solver = client.get_solver()
dwsampler = DWaveSampler(config_file=config_file)
```

Lo más seguro es que nuestro grafo, en cuanto sea un poco grande, no se mapee correctamente con la architectura de qubits del ordenador cúuántico (para valores `length=300, read_length=150, overlap=50`) ya no funciona, asi que buscamos el mapeado usando `minorminer` para obtener el diccionario del grafo embebido:

```
edgelist = solver.edges
adjdict = edgelist_to_adjacency(edgelist)
embed = minorminer.find_embedding(Q, edgelist)
Q_embeded = embed_qubo(Q, embed, adjdict)
```

Resolvemos el problema en el ordenador cuántico y cerramos el cliente:

```
response_qpt = dwsampler.sample_qubo(Q_embeded, num_reads=num_reads)
client.close()
```

Para entender las soluciones tenemos que transformarlas de vuelta desde el mapeado en el grafo embebido:

```
bqm = dimod.BinaryQuadraticModel.from_qubo(Q)
unembedded = unembed_sampleset(response_qpt, embed, bqm, chain_break_method=majority_vote)
```

Ya tenemos las soluciones! Las ordenamos de mayor a menor energía y las reformateamos (de ser un vector de bools volvemos a convertirlas en el formato de `quboDict`):

```
unformated_solutions_list = sorted(embeded_solutions_list.record, key=lambda x: +x[1])
solutions_list = []
for sol, energy, num_appereances in unformated_solutions_list:
	solutions_list = (rebuild_quboDict_from_vector(sol, len(reads)), energy, num_appereances)
```

Como podemos ver aquí las soluciones tienen tres campos: el diccionario con la solución en si, su energía y el número de veces que se ha obtenido esa solución. Entraremos en más detalle en esa frecuencia más adelante. Ejemplo de soluciones obtenidas (aunque aún no puedo comprobar si están bien):

```
Maximum Sampled Configurations from D-Wave	===>
([0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0], -7.45409255, 3)
([0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0], -7.98113883, 2)
([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], -7.98113883, 1)
([0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], -7.45409255, 1)
([0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], -7.45409255, 1)
([0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0], -5.48516016, 1)
([0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], -5.48516016, 1)
([0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0], -5.3797509, 1)
([0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], -5.27434165, 1)
([0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0], -5.16893239, 1)
Minimum Energy Configurations from D-Wave	===>
([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], -7.98113883, 1)
([0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0], -7.98113883, 2)
([0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0], -7.45409255, 3)
([0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], -7.45409255, 1)
([0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0], -7.45409255, 1)
([0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0], -5.48516016, 1)
([0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], -5.48516016, 1)
([0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0], -5.3797509, 1)
([0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], -5.27434165, 1)
([0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0], -5.16893239, 1)
```

Tiempo en ordenador cuántico: Aproximadamente `0.005` segundos.

# Test 2

En el primer test me di cuenta de que la reconstrucción de las solcuiones no era correcta (tomaba siempre las lecturas en el orden lexicográfico: 0, 1, 2, 3... Y como no se barajaban, se obtenía la respuesta correcta).

Una vez subsanado el problema me he dado cuenta de que no siempre tomamos como solución la que empieza por la primera lectura. Recordemos que la solución al QUBO transformado a TSP es un ciclo y no un camino, luego el punto de corte para convertir el ciclo en un camino es determinante. Dos ideas que merece la pena explorar:

- Una vez tenemos el ciclo, cortar por el lado del grafo que tiene mayor energía. Así deberíamos obtener el gen original.

- En vez de buscar un punto de corte, constuir la solución como esté y calcular la [distancia Levenshtein](https://en.wikipedia.org/wiki/Levenshtein_distance#:~:text=Informally%2C%20the%20Levenshtein%20distance%20between,considered%20this%20distance%20in%201965.) (o cualquier otra que nos sirva para calcular distancia entre strings, si existe alguna que compare teniendo en cuenta que principio y final están unidos mejor que mejor), para medir como de cercanas son las solucines.

En el paper de QuASer asumen que las lecturas proceden de un gen circular, asi que quizás este punto no merezca tanto la pena.
