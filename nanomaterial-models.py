# ana.rebeka.kamsek@ki.si, 2023

import numpy as np


def graphene_layers(unit_cells):
    """Constructs graphene layers with dimensions, specified as numbers of unit cells in each direction.

    Graphene layers are randomly rotated and displaced with respect to one another.
    :param unit_cells: numbers of unit cells in x, y, and z
    :return: spatial coordinates of atoms
    """

    # lattice parameters in A
    a = 2.47
    b = 2.47
    c = 3.395

    # how many unit cells
    n_a, n_b, n_c = unit_cells

    # constructing the lattice
    graphene = []

    for k in range(n_c):
        # choose a random rotation angle, less than 10 degrees
        phi = np.random.rand(1) * 10
        phi = np.radians(phi)
        # choose a random displacement, smaller than unit cell parameters
        displacement = (np.random.rand(1) * a, np.random.rand(1) * b)

        for j in range(n_b):
            for i in range(n_a):
                x = i * a * 1 + j * b * 0.5 + k * c * 0
                y = i * a * 0 + j * b * np.sqrt(3) * 0.5 + k * c * 0
                z = i * a * 0 + j * b * 0 + k * c * 1

                x_rotated = float(np.cos(phi) * x + np.sin(phi) * y + displacement[0])
                y_rotated = float(- np.sin(phi) * x + np.cos(phi) * y + displacement[1])
                atom = np.array((x_rotated, y_rotated, z))

                # honeycomb lattice is made up of two displaced hexagonal lattices
                graphene.append(atom)
                graphene.append(atom + (a * 0.5, - b * np.sqrt(3) / 6, 0))

    return np.asarray(graphene)


def create_sphere(element, x, y, z, r=0):
    """Crops a full sphere out of the passed model, assuming a convex model without vacancies.

    :param element: array with element values
    :param x: array with x coordinates of atoms
    :param y: array with y coordinates of atoms
    :param z: array with z coordinates of atoms
    :param r: radius of the desired sphere
    :return: element, x, y, and z values of the sphere
    """

    # in case no radius is supplied, it is determined based on the largest possible full sphere
    # provided the input structure is convex
    if r == 0:
        dimensions = [np.amax(x) - np.amin(x), np.amax(y) - np.amin(y), np.amax(z) - np.amin(z)]
        r = (np.amin(dimensions)) / 2

    # center the structure
    minima = [np.amin(x), np.amin(y), np.amin(z)]
    x = x - minima[0] - r
    y = y - minima[1] - r
    z = z - minima[2] - r

    # crop to yield a sphere
    indices_to_delete = []
    for i in range(len(element)):
        if x[i] ** 2 + y[i] ** 2 + z[i] ** 2 > r ** 2:
            indices_to_delete.append(i)

    element = np.delete(element, indices_to_delete)
    x = np.delete(x, indices_to_delete)
    y = np.delete(y, indices_to_delete)
    z = np.delete(z, indices_to_delete)

    return element, x, y, z


def supported_particle(support, particle, angle=np.pi/4):
    """An option to join two models representing the support and the particle into one model.

    One of the possible configurations for imaging of supported fcc nanoparticles.
    :param support: list of arrays with x, y, and z values for the support
    :param particle: list of arrays with x, y, and z values for the particle
    :param angle: nanoparticle rotation angle, 0 for [100] imaging, np.pi/4 for [110] imaging
    :return: a model of a supported particle
    """

    x_s, y_s, z_s = support
    x_p, y_p, z_p = particle

    # center the support
    x_s = x_s - np.ones(x_s.shape) * (np.amax(x_s) - np.amin(x_s))
    y_s = y_s - np.ones(y_s.shape) * (np.amax(y_s) - np.amin(y_s))
    z_s = -z_s
    support_new = np.vstack((x_s, y_s, z_s)).T

    # rotate the particle
    x_rot = x_p * np.cos(angle) - z_p * np.sin(angle)
    y_rot = y_p
    z_rot = x_p * np.sin(angle) + z_p * np.cos(angle)
    z_rot += np.abs(np.amin(z_rot)) * np.ones(z_rot.shape)
    particle_new = np.column_stack((x_rot, y_rot, z_rot))

    # offset the particle x and y coordinates
    r_particle = np.amax(support_new)
    particle_new[:, 0] += (np.amax(support_new[:, 0]) - r_particle / 2) / 2 * np.ones(particle_new[:, 0].shape)
    particle_new[:, 1] += (np.amax(support_new[:, 1]) - r_particle / 2) / 2 * np.ones(particle_new[:, 1].shape)

    return support_new, particle_new


def write_to_file(filename, element, x, y, z):
    """Creates an .xyz file with atom elements and spatial coordinates.

    :param filename: path to the new file
    :param element: array with element values
    :param x: array with x coordinates of atoms
    :param y: array with y coordinates of atoms
    :param z: array with z coordinates of atoms
    """

    unique_elements = np.unique(element)

    with open(filename, 'w') as new_file:
        new_file.write('{}\n'.format(len(element)))
        for i in range(len(unique_elements)):
            new_file.write('{}'.format(unique_elements[i]))
        new_file.write('\n')

        for i in range(len(element)):
            new_file.write('{} {} {} {}\n'.format(element[i], x[i], y[i], z[i]))
    new_file.close()


def randomize(element, fraction):
    """Assigns random element values among the two in the original element array.

    :param element: array with two distinct element values
    :param fraction: desired fraction of the element that comes first in alphabetical order
    :return: array of randomized element values
    """

    # unique elements are sorted lexicographically, all of them start with a capital letter
    unique_elements = np.unique(element)

    element_randomized = []
    for i in range(len(element)):
        random_choice = np.random.rand(1)
        if random_choice <= fraction:
            element_randomized.append(unique_elements[0])
        else:
            element_randomized.append(unique_elements[1])

    return np.asarray(element_randomized)
