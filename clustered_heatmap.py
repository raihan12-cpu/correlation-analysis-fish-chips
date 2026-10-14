# ============================================================
# USER SETTINGS / EDIT HERE
# ============================================================
from pathlib import Path
from correlation_core import require
ROOT = Path(__file__).resolve().parent
INPUT = ROOT / 'results/analysis_results.json'
OUTPUT = ROOT / 'results/clustered_heatmap'
DPI = 350
ALPHA = 0.05
LABELS = ['Protein', 'Ash', 'Fat', 'Moisture', 'Expansion', 'Oil index', 'Specific volume', 'L*', 'a*', 'b*']
import json
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.lib.colors import Color, HexColor
import pypdfium2 as pdfium

def average_linkage(distance):
    """Agglomerate by mean original pairwise distance; deterministic ties."""
    n = len(distance)
    members = {i: [i] for i in range(n)}
    links = []
    while len(members) > 1:
        ids = sorted(members)
        height, left, right = min(((float(distance[np.ix_(members[a], members[b])].mean()), a, b) for k, a in enumerate(ids) for b in ids[k + 1:]))
        joined = members.pop(left) + members.pop(right)
        links.append([left, right, height, len(joined)])
        members[n + len(links) - 1] = joined
    return (links, next(iter(members.values())))

def plot(matrix, q, links, order, stars):
    """Draw matched row/column dendrograms and annotated correlation cells."""
    name = 'pearson_clustered_BH_stars' if stars else 'pearson_clustered_no_stars'
    path = OUTPUT / (name + '.pdf')
    c = canvas.Canvas(str(path), pagesize=(900, 850))
    c.setTitle('Fish-cracker clustered Pearson correlations')
    x0, y0, cell = (255, 160, 48)
    size = cell * 10
    dark = HexColor('#334655')
    c.setFillColor(dark)
    c.setFont('Helvetica-Bold', 16)
    c.drawString(45, 816, 'Fish-cracker quality attributes')
    c.setFont('Helvetica', 10)
    c.drawString(45, 797, 'Pearson correlations | A-G | 21 independent batches | Carbohydrate excluded')
    pos = {leaf: k for k, leaf in enumerate(order)}
    coords = {i: (pos[i] + 0.5, 0.0) for i in range(10)}
    maxheight = links[-1][2]
    c.setStrokeColor(dark)
    c.setLineWidth(1.5)
    for k, (left, right, height, count) in enumerate(links):
        p1, h1 = coords[left]
        p2, h2 = coords[right]
        coords[10 + k] = ((p1 + p2) / 2, height)
        for orient in ['top', 'left']:

            def point(p, h):
                if orient == 'top':
                    return (x0 + p * cell, y0 + size + 12 + h / maxheight * 125)
                return (140 - h / maxheight * 115, y0 + size - p * cell)
            pts = [point(p1, h1), point(p1, height), point(p2, height), point(p2, h2)]
            for a, b in zip(pts[:-1], pts[1:]):
                c.line(*a, *b)
    for row, i in enumerate(order):
        for col, j in enumerate(order):
            r = float(matrix[i, j])
            value = abs(r)
            end = np.array([0.32, 0.16, 0.52]) if r < 0 else np.array([0.72, 0.34, 0.02])
            rgb = (1 - value) * np.array([0.985, 0.98, 0.975]) + value * end
            c.setFillColor(Color(*rgb))
            c.rect(x0 + col * cell, y0 + (9 - row) * cell, cell, cell, stroke=0, fill=1)
            c.setFillColor(HexColor('#FFFFFF') if value > 0.64 else dark)
            c.setFont('Helvetica', 11)
            txt = f'{(0 if abs(r) < 0.005 else r):.2f}' + ('*' if stars and i != j and (q[i, j] < ALPHA) else '')
            c.drawCentredString(x0 + (col + 0.5) * cell, y0 + (9 - row + 0.5) * cell - 4, txt)
        c.setFillColor(dark)
        c.setFont('Helvetica', 11)
        c.drawRightString(x0 - 10, y0 + (9 - row + 0.5) * cell - 4, LABELS[i])
        c.saveState()
        c.translate(x0 + (row + 0.5) * cell, y0 - 15)
        c.rotate(-45)
        c.drawString(0, 0, LABELS[i])
        c.restoreState()
    bx = 780
    for k in range(500):
        r = -1 + 2 * k / 499
        end = np.array([0.32, 0.16, 0.52]) if r < 0 else np.array([0.72, 0.34, 0.02])
        c.setFillColor(Color(*(1 - abs(r)) * np.array([0.985, 0.98, 0.975]) + abs(r) * end))
        c.rect(bx, y0 + k * size / 500, 17, size / 500 + 0.1, stroke=0, fill=1)
    c.setFillColor(dark)
    c.setFont('Helvetica', 11)
    for r in [-1, -0.5, 0, 0.5, 1]:
        c.drawString(bx + 25, y0 + (r + 1) * size / 2 - 4, f'{r:g}')
    c.saveState()
    c.translate(855, y0 + size / 2)
    c.rotate(90)
    c.drawCentredString(0, 0, 'Pearson r')
    c.restoreState()
    c.setFont('Helvetica', 9)
    c.drawString(45, 50, 'Average linkage; distance = 1 - r. Identical variable order on both axes; branches are exploratory.')
    c.drawString(45, 36, '* BH-adjusted p < 0.05 across 45 unique pairs.' if stars else 'No significance symbols. See numerical results for BH-adjusted p-values.')
    c.drawString(45, 22, 'Oil index: weight-change-based frying measure. Pooled associations can reflect formulation differences.')
    c.save()
    doc = pdfium.PdfDocument(str(path))
    page = doc[0]
    bmp = page.render(scale=DPI / 72)
    bmp.to_pil().save(OUTPUT / (name + '.png'), dpi=(DPI, DPI))
    bmp.close()
    page.close()
    doc.close()

def main():
    report = json.loads(INPUT.read_text())
    matrix = np.array(report['matrices']['pearson'])
    variables = report['audit']['variables']
    q = np.ones((10, 10))
    require(matrix.shape == (10, 10) and np.allclose(matrix, matrix.T), 'Validation failed: matrix.shape == (10, 10) and np.allclose(matrix, matrix.T)')
    require(np.allclose(np.diag(matrix), 1) and np.isfinite(matrix).all(), 'Validation failed: np.allclose(np.diag(matrix), 1) and np.isfinite(matrix).all()')
    for row in report['pearson']:
        i, j = (variables.index(row['x']), variables.index(row['y']))
        q[i, j] = q[j, i] = row['q']
    distance = np.clip(1 - matrix, 0, 2)
    np.fill_diagonal(distance, 0)
    links, order = average_linkage(distance)
    require(sorted(order) == list(range(10)) and len(links) == 9, 'Validation failed: sorted(order) == list(range(10)) and len(links) == 9')
    require(np.all(np.diff([v[2] for v in links]) >= -1e-12), 'Validation failed: np.all(np.diff([v[2] for v in links]) >= -1e-12)')
    test, _ = average_linkage(np.array([[0, 0.2, 0.8], [0.2, 0, 0.6], [0.8, 0.6, 0]]))
    require(np.allclose([r[2] for r in test], [0.2, 0.7]), 'Validation failed: np.allclose([r[2] for r in test], [0.2, 0.7])')
    require(len(report['pearson']) == 45, 'Expected 45 unique pairs')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for stars in [True, False]:
        plot(matrix, q, links, order, stars)
    (OUTPUT / 'clustering_audit.json').write_text(json.dumps({'method': 'average linkage', 'distance': '1 - Pearson r', 'order': [variables[i] for i in order], 'linkage': links, 'significant_unique_pairs': int(np.count_nonzero(q < ALPHA) // 2), 'n': 21, 'known_answer_test': 'passed'}, indent=2))
    print('Saved two PNG/PDF variants; linkage test passed.')
if __name__ == '__main__':
    main()
