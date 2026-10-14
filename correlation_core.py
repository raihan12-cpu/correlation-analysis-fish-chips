"""Validated, import-safe numerical functions for batch-level correlations."""
import itertools
import numpy as np
from scipy.special import betainc


def require(condition, message='Data validation failed'):
    """Raise a clear error even when Python runs with optimization enabled."""
    if not condition:
        raise ValueError(message)


def worksheet_records(sheet, required):
    """Read named columns; reject partial rows and duplicate/missing headers."""
    rows = list(sheet.values)
    require(bool(rows), f'{sheet.title}: empty worksheet')
    headers = list(rows[0])
    for name in required:
        require(headers.count(name) == 1, f'{sheet.title}: expected one column {name}')
    indices = [headers.index(name) for name in required]
    result = []
    for number, row in enumerate(rows[1:], 2):
        values = [row[i] if i < len(row) else None for i in indices]
        if all(value is None for value in values):
            continue  # Fully empty analysis rows or note-only rows.
        require(all(value is not None for value in values), f'{sheet.title}: missing value in row {number}')
        result.append(dict(zip(required, values)))
    require(bool(result), f'{sheet.title}: no observations')
    return result


def keyed(records, fields):
    """Index records without silently overwriting duplicate experimental units."""
    result = {}
    for row in records:
        key = tuple(row[name] for name in fields)
        require(key not in result, f'Duplicate key: {key}')
        result[key] = row
    return result


def ranks(values):
    """Average ranks, including ties."""
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind='stable')
    out = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        out[order[i:j]] = (i + j + 1) / 2
        i = j
    return out


def bh(probabilities):
    """Benjamini-Hochberg adjustment over the supplied comparison family."""
    p = np.asarray(probabilities, dtype=float)
    require(p.ndim == 1 and p.size > 0 and np.isfinite(p).all(), 'Expected finite probability vector')
    require(((p >= 0) & (p <= 1)).all(), 'Probabilities must be within [0, 1]')
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    adjusted[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1,len(p)+1))[::-1])[::-1])
    return adjusted


def correlations(matrix, method, variables):
    """Complete-data correlations; Spearman probabilities are t approximations.

    Missing or constant data fail explicitly. Therefore each pair has the same n.
    BH is applied to unique off-diagonal pairs only.
    """
    data = np.asarray(matrix, dtype=float)
    require(data.ndim == 2 and data.shape[0] >= 3 and data.shape[1] >= 2, 'Need >=3 rows and >=2 variables')
    require(len(variables) == data.shape[1] and len(set(variables)) == len(variables), 'Variable labels do not match unique columns')
    require(np.isfinite(data).all(), 'Missing/nonfinite measurements; no observations silently removed')
    require((np.std(data,axis=0,ddof=1) > 0).all(), 'Constant variable has undefined correlation')
    require(method in ('Pearson', 'Spearman'), 'Unknown correlation method')
    values = np.column_stack([ranks(c) for c in data.T]) if method == 'Spearman' else data
    corr = np.clip(np.corrcoef(values,rowvar=False),-1,1)
    n = len(data)
    result = []
    for i,j in itertools.combinations(range(len(variables)),2):
        coefficient = float(corr[i,j])
        probability = float(betainc((n-2)/2,.5,max(0,1-coefficient*coefficient)))
        result.append(dict(x=variables[i],y=variables[j],r=coefficient,n=n,p=probability))
    for row,adjusted in zip(result,bh([r['p'] for r in result])):
        row['q'] = float(adjusted)
    return corr,result
