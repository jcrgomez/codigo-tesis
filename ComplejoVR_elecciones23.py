# Reimportar librerías necesarias
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from factor_analyzer import FactorAnalyzer
import gudhi as gd
from sklearn.metrics.pairwise import euclidean_distances

# Cargar datos desde el archivo CSV
file_path = "CHES_2024_final_v3.csv"
data = pd.read_csv(file_path, delimiter=";")

# Seleccionar variables relevantes para FA (eliminando variables no significativas)
variables = ['eu_position','lrecon', 'galtan', 'lrgen', 'immigrate_policy',
             'multiculturalism', 'redistribution', 'climate_change',
             'environment', 'nationalism', 'protectionism', 'regions',
             'deregulation', 'womens_rights', 'lgbtq_rights', 'samesex_marriage',
             'religious_principles', 'ethnic_minorities', 'executive_power', 'judicial_independence',
             'corrupt_salience', 'people_v_elite']

# Limpiar datos y estandarizar
data_clean = data[['party'] + variables].dropna()
scaler = StandardScaler()
scaled_data = scaler.fit_transform(data_clean[variables])

# Aplicar Análisis Factorial (FA)
fa = FactorAnalyzer(n_factors=2, rotation="varimax")
fa.fit(scaled_data)

# Obtener cargas factoriales
factor_loadings = pd.DataFrame(fa.loadings_, index=variables, columns=['Factor 1', 'Factor 2'])

# Transformar los datos originales usando FA
fa_components = fa.transform(scaled_data)
fa_df = pd.DataFrame(fa_components, columns=['Factor 1', 'Factor 2'])
fa_df['party'] = data_clean['party'].values

# Calcular la matriz de distancias en el espacio de los factores
coords = fa_df[['Factor 1', 'Factor 2']].values
distance_matrix = euclidean_distances(coords, coords)

# Modificar la dimensión máxima de los símplices
max_simplex_size = 8  # Cambia este valor para permitir símplices de mayor tamaño

# Construcción del complejo de Vietoris-Rips con GUDHI
epsilon = 2.1  # Umbral de conexión
rips_complex = gd.RipsComplex(points=coords, max_edge_length=epsilon)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_simplex_size)

# Extraer los símplices de la dimensión deseada
simplices = [simplex[0] for simplex in simplex_tree.get_skeleton(max_simplex_size)]
simplices = [s for s in simplices if len(s) <= max_simplex_size]  # Filtrar símplices grandes

# Convertir índices a nombres de partidos
party_names = fa_df['party'].values
complejo = {frozenset([party_names[i] for i in simplex]) for simplex in simplices}

# Filtrar los símplices que contienen a PP y PSOE simultáneamente
complejo = {s for s in complejo if not ("PP" in s and "PSOE" in s)}

complejo = {s for s in complejo if not ("Vox" in s and "Junts" in s)}

# 📊 Graficar el Análisis Factorial (FA) con SOLO la envoltura convexa de los símplices
plt.figure(figsize=(8, 6))

# Dibujar puntos con nombres de partidos
for _, row in fa_df.iterrows():
    plt.scatter(row['Factor 1'], row['Factor 2'], alpha=0.8)
    plt.text(row['Factor 1'] + 0.02, row['Factor 2'] + 0.02, row['party'], fontsize=11, fontweight='bold')

for simplex in complejo:
    vertices = [fa_df[fa_df['party'] == party][['Factor 1', 'Factor 2']].values[0] for party in simplex]
    
    # Dibujar líneas si el simplex es una arista (2 vértices)
    if len(vertices) == 2:
        plt.plot([vertices[0][0], vertices[1][0]], [vertices[0][1], vertices[1][1]], 'k-', linewidth=0.6)

    # Dibujar envoltura convexa si el simplex tiene más de 2 vértices
    if len(vertices) > 2:
        polygon = plt.Polygon(vertices, color='gray', alpha=0.2, edgecolor='k')
        plt.gca().add_patch(polygon)


# Configuración del gráfico
plt.axhline(0, color='gray', linestyle='--', alpha=0.5)
plt.axvline(0, color='gray', linestyle='--', alpha=0.5)
plt.xlabel("Izquierda vs. Derecha Económica",fontsize=11)
plt.ylabel("Institucionalismo vs.Populismo / UE-Soberanismo",fontsize=11)
plt.title(f" Complejo Vietoris-Rips (ε = {epsilon})",fontsize=15)

plt.show()

