from pathlib import Path

import numpy as np

from hsqc_methods import nn_peak_similarity, quadtree_similarity
from hsqc_similarity import Spectrum2D


def gaussian_2d(f2, f1, c2, c1, w2, w1, height=1.0):
    g2 = np.exp(-0.5 * ((f2 - c2) / w2) ** 2)
    g1 = np.exp(-0.5 * ((f1 - c1) / w1) ** 2)
    return height * np.outer(g1, g2)


def make_spectrum(peaks, source):
    ppm_f2 = np.linspace(0, 10, 256)
    ppm_f1 = np.linspace(0, 160, 256)
    intensity = np.zeros((ppm_f1.size, ppm_f2.size))
    for c2, c1, w2, w1, h in peaks:
        intensity += gaussian_2d(ppm_f2, ppm_f1, c2, c1, w2, w1, h)
    return Spectrum2D(ppm_f2=ppm_f2, ppm_f1=ppm_f1, intensity=intensity, source=Path(source))


REF = [(3.0, 40.0, 0.1, 1.5, 1.0), (7.0, 120.0, 0.1, 1.5, 0.8)]
SHIFTED = [(3.1, 41.0, 0.1, 1.5, 1.0), (7.1, 121.0, 0.1, 1.5, 0.8)]
UNRELATED = [(5.0, 80.0, 0.1, 1.5, 1.0), (1.0, 20.0, 0.1, 1.5, 0.8)]


def test_quadtree_identity_is_one():
    spectrum = make_spectrum(REF, "synthetic")
    assert quadtree_similarity(spectrum, spectrum)["similarity"] == 1.0


def test_quadtree_related_beats_unrelated():
    reference = make_spectrum(REF, "ref")
    shifted = make_spectrum(SHIFTED, "shifted")
    unrelated = make_spectrum(UNRELATED, "unrelated")
    assert (
        quadtree_similarity(reference, shifted)["similarity"]
        > quadtree_similarity(reference, unrelated)["similarity"]
    )


def test_nn_identity_is_one():
    spectrum = make_spectrum(REF, "synthetic")
    result = nn_peak_similarity(spectrum, spectrum)
    assert result["similarity"] == 1.0
    assert result["mean_distance"] == 0.0


def test_nn_related_beats_unrelated():
    reference = make_spectrum(REF, "ref")
    shifted = make_spectrum(SHIFTED, "shifted")
    unrelated = make_spectrum(UNRELATED, "unrelated")
    assert (
        nn_peak_similarity(reference, shifted)["similarity"]
        > nn_peak_similarity(reference, unrelated)["similarity"]
    )
