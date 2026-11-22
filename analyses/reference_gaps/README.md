# Análisis in silico de puntos ciegos del eDNA para criptobentónicos del Golfo de California

Proyecto CBF-2026-4124. Análisis preliminar del 21 de noviembre de 2026. Los datos de GenBank y OBIS cambian con el tiempo; el script `cripticos_insilico.py` reproduce el análisis con la fecha de consulta en que se ejecute.

## Archivos

| Archivo | Contenido |
| --- | --- |
| `cripticos_insilico.py` | Script reproducible (Python 3, `requests`) |
| `results/` (no incluido) | Tablas de resultados; se publicarán con el artículo correspondiente |

## Métodos

**Universo de especies.**
- Lista de especies desde la API de OBIS (`/v3/checklist`), con un polígono del Golfo de California. El polígono excluye la costa pacífica de la península.
- Solo se incluyeron especies con al menos 3 registros.
- Criptobentónicos: las 9 familias del protocolo (Apogonidae, Blenniidae, Bythitidae, Chaenopsidae, Dactyloscopidae, Gobiesocidae, Gobiidae, Labrisomidae, Tripterygiidae).
- Comparador conspicuo: 15 familias arrecifales (Pomacentridae, Labridae, Scaridae, Serranidae, Epinephelidae, Lutjanidae, Haemulidae, Acanthuridae, Chaetodontidae, Pomacanthidae, Balistidae, Tetraodontidae, Mullidae, Kyphosidae, Cirrhitidae).
- Se añadieron 5 especies criptobentónicas de la tesis de El Khattabi Salazar (2014) sin registros en OBIS dentro del polígono.

**Referencias.**
- Búsqueda en GenBank (nuccore) por especie, a través de la taxonomía de NCBI, de modo que se incluyen sinónimos.
- Consultas: 12S o mitogenomas (hasta 30 registros por especie), COI (conteo) y mitogenomas completos (conteo).
- Cada secuencia se ubicó en coordenadas del 12S del mitogenoma de *Holacanthus passer* (NC_027594.1) mediante k-meros únicos (k = 8) y votación diagonal.
- Una región de marcador (el inserto entre cebadores) se consideró:
  - "completa" si la secuencia la abarca con una tolerancia de 20 pb en los extremos;
  - "parcial" si cubre al menos el 50 % del inserto.

**Cebadores.**
- MiFish-U (Miya et al., 2015), teleo (Valentini et al., 2016) y Riaz 12S-V5 (Riaz et al., 2011).
- Se contó el número de desajustes en el mejor sitio, sin huecos y en ambas hebras, y cuántos caen en las 5 bases del extremo 3′.
- Por especie solo se usaron secuencias que contienen el sitio: mejor sitio con 6 desajustes o menos.
- Por familia se usaron hasta 60 mitogenomas completos de GenBank de cualquier región.

**Resolución.**
- Se extrajo el inserto entre cebadores cuando ambos sitios tienen 6 desajustes o menos.
- Para cada especie se calculó la distancia de edición mínima a cualquier otra especie del universo con amplicón disponible.
- Distancia 0 = amplicón idéntico a otra especie (no identificable). Distancia 1 a 2 = identificación frágil ante errores de secuenciación o variación intraespecífica.

## Limitaciones

- **Lista de especies:** OBIS contiene registros mal identificados o fuera de rango; el umbral de 3 registros los reduce pero no los elimina.
- **Taxonomía:** la familia viene de la taxonomía de OBIS, que puede diferir de la del protocolo (por ejemplo, Epinephelidae separado de Serranidae).
- **Sinónimos:** NCBI puede no reconocer sinónimos recientes, lo que subestima la cobertura de algunas especies.
- **Sesgo de los mitogenomas:** el análisis por familia usa mitogenomas de cualquier región y no necesariamente especies del Golfo.
- **Desajustes y PCR:** el número de desajustes no predice por sí solo la eficiencia real de la PCR, que depende de su posición, del tipo de desajuste y de las condiciones de reacción.
- **Coordenadas:** la ubicación por k-meros es aproximada para secuencias cortas o divergentes.
- **Errores en GenBank:** amplicones idénticos entre especies pueden reflejar errores de identificación en GenBank.

## Mitogenomes from public raw reads
See [`sra_mitogenomes/`](sra_mitogenomes/): recovery of mitochondrial references from SRA/ENA raw reads (WGS and targeted capture) to fill marker gaps for Gulf of California fishes.
