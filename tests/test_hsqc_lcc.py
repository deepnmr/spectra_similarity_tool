from pathlib import Path

import numpy as np

from hsqc_lcc import lcc_similarity
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

KW = dict(sigma_f2=0.1, sigma_f1=1.5, step_f2=0.04, step_f1=0.6)


def test_lcc_identity_is_one():
    spectrum = make_spectrum(REF, "synthetic")
    assert lcc_similarity(spectrum, spectrum, **KW)["similarity"] == 1.0


def test_lcc_related_beats_unrelated():
    reference = make_spectrum(REF, "ref")
    shifted = make_spectrum(SHIFTED, "shifted")
    unrelated = make_spectrum(UNRELATED, "unrelated")
    assert (
        lcc_similarity(reference, shifted, **KW)["similarity"]
        > lcc_similarity(reference, unrelated, **KW)["similarity"]
    )


def test_lcc_monotonic_in_drift():
    # Graded shift tolerance: growing drift lowers the score monotonically, no bin-edge jumps.
    reference = make_spectrum(REF, "ref")
    scores = []
    for d in (0.0, 0.05, 0.1, 0.2, 0.4):
        drifted = make_spectrum([(3.0 + d, 40.0 + 10 * d, 0.1, 1.5, 1.0),
                                 (7.0 + d, 120.0 + 10 * d, 0.1, 1.5, 0.8)], "drift")
        scores.append(lcc_similarity(reference, drifted, **KW)["similarity"])
    assert all(scores[i] >= scores[i + 1] - 1e-9 for i in range(len(scores) - 1)), scores
