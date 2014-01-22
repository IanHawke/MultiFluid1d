import numpy as np

def EOS_Paper(n, options):
    """
    This encodes the master function used in the original two-stream paper.

    Parameters
    ----------

    n : float, array
        A 2 component [t, x], 2 species [1, 2] array with the background for the number current.
    options : dict
              All other options        
    """

    # n_{XY}^2 = -g_{mu nu} n_{X}^{mu} n_{Y}^mu

    nsq = np.zeros(3)
    nsq[0] = n[0, 0]**2 - n[0, 1]**2
    nsq[1] = n[1, 0]**2 - n[1, 1]**2
    nsq[2] = n[0, 0] * n[1, 0] - n[0, 1] * n[1, 1]
    nn = np.sqrt(nsq)

    dLdn = np.zeros_like(nn)
    d2Ldn = np.zeros((len(nn), len(nn)))
    
    m = options['m']
    kappa = options['kappa']
    kappa_12 = options['kappa_12']
    kappa_d = options['kappa_delta']
    sigma = options['sigma']
    g = options['gamma']
    
#    L = -m[0] * nn[0] - m[1] * nn[1] - kappa[0] * nn[0] ** g[0] - kappa[1] * nn[1] ** g[1] - kappa_12 * nn[0] ** sigma[0] * nn[1] ** sigma[1] - kappa_d * (1 - nn[0] ** 2 * nn[1] ** 2 / nn[2] ** 2)
    dLdn[0] = -(1 / nn[2] ** 2 / nn[0] * (m[0] * nn[2] ** 2 + kappa[0] * nn[0] ** (g[0] - 1) * g[0] * nn[2] ** 2 + kappa_12 * nn[0] ** (sigma[0] - 1) * sigma[0] * nn[1] ** sigma[1] * nn[2] ** 2 - 2 * kappa_d * nn[0] * nn[1] ** 2)) / 2.0
    dLdn[1] = -(1 / nn[2] ** 2 / nn[1] * (m[1] * nn[2] ** 2 + kappa[1] * nn[1] ** (g[1] - 1) * g[1] * nn[2] ** 2 + kappa_12 * nn[0] ** sigma[0] * nn[1] ** (sigma[1] - 1) * sigma[1] * nn[2] ** 2 - 2 * kappa_d * nn[0] ** 2 * nn[1])) / 2.0
    dLdn[2] = -kappa_d * nn[0] ** 2 * nn[1] ** 2 / nn[2] ** 4
    d2Ldn[0, 0] = -(1 / nn[0] ** 3 * (nn[0] ** (g[0] - 1) * kappa[0] * g[0] ** 2 - 2 * kappa[0] * nn[0] ** (g[0] - 1) * g[0] + nn[0] ** (sigma[0] - 1) * kappa_12 * sigma[0] ** 2 * nn[1] ** sigma[1] - 2 * kappa_12 * nn[0] ** (sigma[0] - 1) * sigma[0] * nn[1] ** sigma[1] - m[0])) / 4.0
    d2Ldn[1, 1] = -(1 / nn[1] ** 3 * (nn[1] ** (g[1] - 1) * kappa[1] * g[1] ** 2 - 2 * kappa[1] * nn[1] ** (g[1] - 1) * g[1] + nn[1] ** (sigma[1] - 1) * kappa_12 * nn[0] ** sigma[0] * sigma[1] ** 2 - 2 * kappa_12 * nn[0] ** sigma[0] * nn[1] ** (sigma[1] - 1) * sigma[1] - m[1])) / 4.0
    d2Ldn[2, 2] = 2 * kappa_d * nn[0] ** 2 * nn[1] ** 2 / nn[2] ** 6
    d2Ldn[0, 1] = -(1 / nn[2] ** 2 / nn[0] / nn[1] * (kappa_12 * nn[0] ** (sigma[0] - 1) * sigma[0] * nn[1] ** (sigma[1] - 1) * sigma[1] * nn[2] ** 2 - 4 * kappa_d * nn[0] * nn[1])) / 4.0
    d2Ldn[0, 2] = -kappa_d * nn[1] ** 2 / nn[2] ** 4
    d2Ldn[1, 2] = -kappa_d * nn[0] ** 2 / nn[2] ** 4
    d2Ldn[1, 0] = d2Ldn[0, 1]
    d2Ldn[2, 0] = d2Ldn[0, 2]
    d2Ldn[2, 1] = d2Ldn[1, 2]

    return (dLdn, d2Ldn)
    
