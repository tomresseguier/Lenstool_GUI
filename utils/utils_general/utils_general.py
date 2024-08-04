import numpy as np
from scipy.spatial.distance import pdist, squareform



def read_text_file(file_path):
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            cleaned_line = line.strip()
            if not cleaned_line.startswith("#"):
                lines.append(cleaned_line)
    return lines


def find_close_coord(cat, d):
    ids = cat['id']
    x = cat['x']
    y = cat['y']
    
    points = np.vstack((x, y)).T
    distance_matrix = squareform(pdist(points))
    n = len(points)
    
    parent = list(range(n))
    
    def find(u):
        if parent[u] != u:
            parent[u] = find(parent[u])
        return parent[u]
    
    def union(u, v):
        root_u = find(u)
        root_v = find(v)
        if root_u != root_v:
            parent[root_u] = root_v

    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i, j] < d:
                union(i, j)
    
    groups = {}
    for i in range(n):
        root = find(i)
        if root in groups:
            groups[root].append(ids[i])
        else:
            groups[root] = [ids[i]]
    
    result = list(groups.values())
    
    return result




