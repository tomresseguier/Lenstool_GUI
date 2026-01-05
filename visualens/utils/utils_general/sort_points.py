import numpy as np
from scipy.spatial import cKDTree
from tqdm import tqdm

def unit(v):
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    return v / n if n > 0 else np.zeros_like(v)

def local_tangent_vectors(points, k=12, tree=None):
    """
    Return a (N,2) array of unit tangent direction estimates for each point.
    Uses PCA on k nearest neighbors (including the point).
    """
    pts = np.asarray(points)
    N = len(pts)
    if tree is None:
        tree = cKDTree(pts)
    dists, idxs = tree.query(pts, k=k)
    tangents = np.zeros((N, pts.shape[1]))
    for i in range(N):
        nbrs = pts[idxs[i]]
        # center and do 2D PCA (covariance)
        X = nbrs - nbrs.mean(axis=0)
        cov = X.T @ X
        # eigenvector with largest eigenvalue
        w, v = np.linalg.eigh(cov)
        principal = v[:, np.argmax(w)]
        tangents[i] = unit(principal)
    return tangents

def angle_penalty(prev_dir, new_vec):
    """
    angle penalty term in [0, 1], where 0 = aligned, 1 = opposite/perp extreme.
    We'll use 1 - cos(theta) in [0,2]; scale appropriately in cost.
    """
    u_prev = unit(prev_dir)
    u_new = unit(new_vec)
    cosang = np.clip(np.dot(u_prev, u_new), -1.0, 1.0)
    return 1.0 - cosang  # 0 when aligned, 2 when opposite

def remove_doubles_and_nan(coords) :
    coords_x = np.unique(coords[0])
    coords_y = np.unique(coords[1])
    coords_unique = [coords_x, coords_y] if len(coords_x) == len(coords_y) else coords
    nan_mask_x = ~np.isnan(coords_unique[0])
    nan_mask_y = ~np.isnan(coords_unique[1])
    nan_mask = nan_mask_x & nan_mask_y
    pts = np.asarray(coords_unique, dtype=float).T[nan_mask].T
    return pts

def sort_points(coords,
                k_tangent=12,
                k_candidates=10,
                alpha=1.0,
                distance_threshold=None,
                angle_threshold=np.pi/2):
    """
    coords: (x, y) coordinates of points
    k_tangent: neighbors for local tangent estimation
    k_candidates: how many nearest neighbor candidates to consider when extending
    alpha: multiplier for angular penalty (higher -> prefers direction over short jump)
    distance_threshold: maximum allowed step distance (if None, set to mean+3*std)
    angle_threshold: maximum allowed absolute angle (radians) between prev direction
                     and candidate segment (defaults to 90 degrees)
    Returns: (M,2) array of ordered points with np.nan rows separating lines.
    """
    pts = remove_doubles_and_nan(coords).T
    N = len(pts)
    if N == 0:
        return np.empty((0, 2))
    tree = cKDTree(pts)
    if distance_threshold is None:
        # heuristic: set threshold to something reasonable
        dists, _ = tree.query(pts, k=5)
        # dists includes self=0; take first non-zero neighbor
        nearest_nonzero = dists[:, 1]
        distance_threshold = np.mean(nearest_nonzero) + 3*np.std(nearest_nonzero)

    tangents = local_tangent_vectors(pts, k=k_tangent, tree=tree)

    visited = np.zeros(N, dtype=bool)
    indices_remaining = set(range(N))
    ordered_list = []

    def pick_seed():
        # Prefer points with small number of neighbors within small radius -> endpoints
        # compute neighbor counts within small radius (e.g. median dist)
        med = np.median(tree.query(pts, k=6)[0][:,1])
        counts = tree.query_ball_point(pts, r=med)
        # convert to counts
        counts_len = np.array([len(c) for c in counts])
        # choose unvisited point with minimal neighbor count (endpoint-like)
        cand = np.where(~visited)[0]
        if cand.size == 0:
            return None
        idx = cand[np.argmin(counts_len[cand])]
        return idx

    while len(indices_remaining) > 0:
        seed = pick_seed()
        if seed is None:
            break
        # start a path from seed; we'll grow forwards and backwards
        for direction in [1, -1]:  # first extend forward, then backward
            if visited[seed]:
                break
            path = [seed] if direction == 1 else []
            # If backward: we'll build then reverse and merge; simpler approach: always grow forward from seed,
            # then try to grow in opposite direction starting from seed again to prepend.
            # Simpler: do only forward growth here; after finishing, try to attach at both ends.
            # For clarity we implement grow function and then attach both ends later.
            break

        # Implement grow from a start index and optionally prepend
        def grow(start_idx, prepend=False):
            cur_idx = start_idx
            prev_dir = None
            if prepend and ordered_path:
                prev_dir = unit(pts[ordered_path[0]] - pts[cur_idx])  # initial direction when prepending
            local_path = []
            while True:
                if visited[cur_idx]:
                    # if already visited skip
                    break
                visited[cur_idx] = True
                indices_remaining.discard(cur_idx)
                local_path.append(cur_idx)

                # find nearest unvisited neighbors
                dists, idxs = tree.query(pts[cur_idx], k=min(k_candidates, N))
                # idxs may include already visited or self
                best = None
                best_cost = np.inf
                for d, j in zip(np.atleast_1d(dists), np.atleast_1d(idxs)):
                    if j == cur_idx or visited[j]:
                        continue
                    if d > distance_threshold:
                        continue
                    candidate_vec = pts[j] - pts[cur_idx]
                    # define a preferred direction to compare: if we have prev_dir use it, else use local tangent
                    if prev_dir is None:
                        pref_dir = tangents[cur_idx]
                    else:
                        pref_dir = prev_dir
                    penalty = angle_penalty(pref_dir, candidate_vec)
                    ang = np.arccos(np.clip(np.dot(unit(pref_dir), unit(candidate_vec)), -1, 1))
                    if ang > angle_threshold:
                        continue
                    cost = d + alpha * penalty
                    if cost < best_cost:
                        best_cost = cost
                        best = j
                        best_ang = ang
                        best_vec = candidate_vec
                if best is None:
                    break
                # accept best
                prev_dir = unit(best_vec)
                cur_idx = best
            return local_path

        # Start by growing forward from seed
        ordered_path = grow(seed, prepend=False)
        # Try to grow backwards from seed and prepend (walk from seed but choose previous direction)
        # To grow backwards, we want to consider candidates that lead toward seed: reverse sense by flipping vector sign in angle check.
        # Implement a small wrapper that uses negative prev_dir to bias opposite direction.
        # Simpler: run grow from seed but when selecting candidates, prefer ones that make prev_dir opposite.
        # We'll just run grow from seed with prev_dir = -tangent to encourage opposite direction and then prepend reversed.
        # Set prev_dir initial to -tangents[seed] and modify grow to accept initial prev_dir.
        def grow_with_initial(start_idx, initial_prev_dir):
            cur_idx = start_idx
            prev_dir = unit(initial_prev_dir)
            local = []
            while True:
                if visited[cur_idx]:
                    break
                visited[cur_idx] = True
                indices_remaining.discard(cur_idx)
                local.append(cur_idx)
                dists, idxs = tree.query(pts[cur_idx], k=min(k_candidates, N))
                best = None
                best_cost = np.inf
                for d, j in zip(np.atleast_1d(dists), np.atleast_1d(idxs)):
                    if j == cur_idx or visited[j]:
                        continue
                    if d > distance_threshold:
                        continue
                    candidate_vec = pts[j] - pts[cur_idx]
                    penalty = angle_penalty(prev_dir, candidate_vec)
                    ang = np.arccos(np.clip(np.dot(unit(prev_dir), unit(candidate_vec)), -1, 1))
                    if ang > angle_threshold:
                        continue
                    cost = d + alpha * penalty
                    if cost < best_cost:
                        best_cost = cost
                        best = j
                        best_vec = candidate_vec
                if best is None:
                    break
                prev_dir = unit(best_vec)
                cur_idx = best
            return local

        # We already consumed the seed in ordered_path (visited). To get the other side, unvisit the seed temporarily
        # (so the second grow can include it as starting point) — but we need to be careful with visited flags.
        # Simpler approach: create a second path by running grow_with_initial but without marking visited so we can prepend.
        # We'll run a "lookahead" growth that ignores visited but only returns indices not yet in ordered_path.
        def extend_ignoring_visited(start_idx, initial_prev_dir):
            cur_idx = start_idx
            prev_dir = unit(initial_prev_dir)
            local = []
            seen_local = set()
            while True:
                # don't mark visited here; we will filter later
                local.append(cur_idx)
                seen_local.add(cur_idx)
                dists, idxs = tree.query(pts[cur_idx], k=min(k_candidates, N))
                best = None
                best_cost = np.inf
                for d, j in zip(np.atleast_1d(dists), np.atleast_1d(idxs)):
                    if j == cur_idx or j in seen_local or j in ordered_path:
                        continue
                    if d > distance_threshold:
                        continue
                    candidate_vec = pts[j] - pts[cur_idx]
                    penalty = angle_penalty(prev_dir, candidate_vec)
                    ang = np.arccos(np.clip(np.dot(unit(prev_dir), unit(candidate_vec)), -1, 1))
                    if ang > angle_threshold:
                        continue
                    cost = d + alpha * penalty
                    if cost < best_cost:
                        best_cost = cost
                        best = j
                        best_vec = candidate_vec
                if best is None:
                    break
                prev_dir = unit(best_vec)
                cur_idx = best
            return local  # includes start_idx first

        # already consumed ordered_path (visited flagged). Build backward extension without touching visited flags:
        back_seed_dir = -tangents[seed]  # encourage opposite direction
        backward_guess = extend_ignoring_visited(seed, back_seed_dir)
        # backward_guess starts at seed; discard the seed (will be duplicated) and reverse
        backward_chain = [idx for idx in backward_guess if idx not in ordered_path]
        backward_chain = backward_chain[::-1]  # reverse to prepend
        # now mark backward_chain visited properly and prepend to ordered_path
        for idx in backward_chain:
            if not visited[idx]:
                visited[idx] = True
                indices_remaining.discard(idx)
        full_chain = backward_chain + ordered_path

        # append the chain in coordinates with nan separator
        if len(full_chain) > 0:
            ordered_list.append(pts[full_chain])
        
        print( str( (100 - len(indices_remaining)/N * 100) ) + '% of points sorted' )

    # concatenate with nan rows between lines
    if len(ordered_list) == 0:
        return np.empty((0,2))
    out = []
    for i, arr in enumerate(ordered_list):
        out.append(arr)
        if i != len(ordered_list)-1:
            out.append(np.array([[np.nan, np.nan]]))
    return np.vstack(out).T







def break_curves(coords, distance_threshold=8.) :
    coords_unique = coords#remove_doubles_and_nan(coords)
    breaks = []
    for i in tqdm(range( len(coords_unique[0])-1 )) :
        xi, yi = coords_unique[0][i], coords_unique[1][i]
        xf, yf = coords_unique[0][i+1], coords_unique[1][i+1]
        d = ((xf-xi)**2 + (yf-yi)**2)**0.5
        if d>distance_threshold :
            breaks.append(i+1)
    breaks.reverse()
    for i in breaks :
        coords_unique[0] = np.insert(coords_unique[0], i, np.nan)
        coords_unique[1] = np.insert(coords_unique[1], i, np.nan)
    return coords_unique



