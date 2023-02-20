from sys import argv

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#used if no GUI is available
matplotlib.use('Agg')

#fix figure size
plt.rcParams['figure.figsize'] = [20, 20]

if len(argv) < 2:
    print('usage: python plot_mesh.py mesh.csv outputfilename')
    exit(-1)

vertices = []
segments = []
triangles = []
quads = []
polygons = []

#switch to show entities' indices in plot
annotate = True

with open(argv[1]) as f:
    for line in f:
        array = np.fromstring(line.split()[0], sep=',')

        if array[0] == 0:
            if array.size == 5:
                triangles.append(array)
            elif array.size == 6:
                quads.append(array)
            else:
                polygons.append(array)
        if array[0] == 1:
            segments.append(array)
        if array[0] == 2:
            vertices.append(array)

vertices = np.vstack(vertices)
segments = np.vstack(segments)
if triangles != []:
    triangles = np.vstack(triangles)
if quads != []:
    quads = np.vstack(quads)

# plot vertices
for i, num in enumerate(vertices[:, 1]):
    if(annotate):
            plt.annotate(
                int(num),
                vertices[np.where(vertices[:, 1] == num)][0, -2:],
                fontsize='xx-large',
                weight='bold',
                bbox=dict(boxstyle='circle', facecolor='none', edgecolor='red'),
                color='r',
                ha='center',
                va='center'
            )   

# plot segments
for segment in segments:
    x_coords, y_coords = np.column_stack((
        vertices[np.where(vertices[:, 1] == int(segment[2]))][0, -2:],
        vertices[np.where(vertices[:, 1] == int(segment[3]))][0, -2:]
    ))
    plt.plot(x_coords, y_coords, 'k-')

    midpoint = .5 * (
            vertices[np.where(vertices[:, 1] == segment[2])][0, -2:] +
            vertices[np.where(vertices[:, 1] == segment[3])][0, -2:]
    )
    if(annotate):
        plt.annotate(
            int(segment[1]),
            midpoint,
            fontsize='xx-large',
            ha='center',
            va='center'
        )


# plot cells
def plot_cells(cells):
    for cell in cells:
        cell_vertices = []
        for i in range(2, cell.size):
            segment = segments[np.where(segments[:, 1] == cell[i])][0]
            cell_vertices.append(
                vertices[np.where(vertices[:, 1] == segment[2])][0, -2:]
            )
            cell_vertices.append(
                vertices[np.where(vertices[:, 1] == segment[3])][0, -2:]
            )

        cell_vertices = np.vstack(cell_vertices)
        if(annotate):
            plt.annotate(
                int(cell[1]),
                np.mean(cell_vertices, axis=0),
                fontsize='xx-large',
                color='deeppink',
                weight='bold',
                ha='center',
                va='center'
            )


plot_cells(triangles)
plot_cells(quads)
plot_cells(polygons)

plt.axis('off')

out_file = f"{argv[2]}.png"

#plt.show()
plt.savefig(out_file, bbox_inches='tight')
