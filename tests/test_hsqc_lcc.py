from pathlib import Path

import numpy as np
import pytest

from hsqc_lcc import build_parser, cosine_similarity, lcc_similarity, local_contrast_similarity
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


def test_cosine_self_is_one():
    spectrum = make_spectrum(REF, "synthetic")
    assert cosine_similarity(spectrum, spectrum, **KW)["similarity"] == 1.0


def test_mean_centring_is_the_discriminating_step():
    # The ablation claim: un-centred cosine of all-positive images cannot penalise
    # non-co-located intensity, so it scores an UNRELATED spectrum HIGHER than the
    # mean-centred LCC does. Mean-centring is what turns overlap into discrimination.
    reference = make_spectrum(REF, "ref")
    unrelated = make_spectrum(UNRELATED, "unrelated")
    lcc = lcc_similarity(reference, unrelated, **KW)["similarity"]
    cos = cosine_similarity(reference, unrelated, **KW)["similarity"]
    assert cos > lcc


def test_abs_baseline_keeps_negative_ch2_peaks():
    # Edited-HSQC CH2 peaks are negative; clip deletes them (so a spectrum with an extra
    # negative peak looks identical to one without), abs keeps them (so it looks different).
    a = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0)], "ch")               # one positive peak
    b = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0),
                       (7.0, 120.0, 0.1, 1.5, -0.8)], "ch+ch2")         # + a negative CH2 peak
    clipped = lcc_similarity(a, b, baseline="clip", **KW)["similarity"]
    absed = lcc_similarity(a, b, baseline="abs", **KW)["similarity"]
    assert clipped > absed  # under clip the extra CH2 vanishes, so b looks like a


def test_local_contrast_contract_and_metadata():
    reference = make_spectrum(REF, "ref")
    shifted = make_spectrum(SHIFTED, "shifted")
    unrelated = make_spectrum(UNRELATED, "unrelated")
    scaled = Spectrum2D(
        shifted.ppm_f2, shifted.ppm_f1, shifted.intensity * 17.0, Path("scaled")
    )

    result = local_contrast_similarity(reference, shifted, **KW)
    reverse = local_contrast_similarity(shifted, reference, **KW)["similarity"]
    scores = [
        local_contrast_similarity(reference, reference, **KW)["similarity"],
        result["similarity"], reverse,
        local_contrast_similarity(reference, unrelated, **KW)["similarity"],
    ]
    assert scores[0] == 1.0
    assert all(0.0 <= score <= 1.0 for score in scores)
    assert result["similarity"] == pytest.approx(reverse, abs=1e-12)
    assert result["similarity"] == pytest.approx(
        local_contrast_similarity(reference, scaled, **KW)["similarity"], abs=1e-12
    )
    assert result["method"] == "local-contrast"
    assert result["background_factor"] == 3
    assert result["intensity_transform"] == "sqrt"
    assert result["sigma_f2"] == KW["sigma_f2"]
    assert result["sigma_f1"] == KW["sigma_f1"]
    assert result["range_f2"] == [0.0, 10.0]
    assert result["range_f1"] == [0.0, 160.0]


def test_local_contrast_decreases_monotonically_with_shift():
    reference = make_spectrum(REF, "ref")
    scores = []
    for drift in (0.0, 0.05, 0.1, 0.2, 0.4):
        shifted = make_spectrum([
            (3.0 + drift, 40.0 + 10 * drift, 0.1, 1.5, 1.0),
            (7.0 + drift, 120.0 + 10 * drift, 0.1, 1.5, 0.8),
        ], "shifted")
        scores.append(local_contrast_similarity(reference, shifted, **KW)["similarity"])
    assert all(a >= b - 1e-9 for a, b in zip(scores, scores[1:])), scores


def test_local_contrast_rejects_baseline_and_noise_before_unrelated_spectrum():
    reference = make_spectrum(REF, "ref")
    related = make_spectrum(SHIFTED, "related")
    broad = gaussian_2d(related.ppm_f2, related.ppm_f1, 5.0, 80.0, 4.0, 60.0, 0.03)
    noisy_related = Spectrum2D(
        related.ppm_f2,
        related.ppm_f1,
        related.intensity + broad + np.random.default_rng(0).normal(0, 0.01, related.intensity.shape),
        related.source,
    )
    unrelated = make_spectrum(UNRELATED, "unrelated")
    assert (
        local_contrast_similarity(reference, noisy_related, **KW)["similarity"]
        > local_contrast_similarity(reference, unrelated, **KW)["similarity"]
    )


def test_local_contrast_is_grid_invariant_at_converged_resolution():
    ppm_f2 = np.arange(0.0, 2.0, 0.01)
    ppm_f1 = np.arange(0.0, 20.0, 0.1)

    def fine_spectrum(drift):
        image = sum(
            gaussian_2d(ppm_f2, ppm_f1, c2 + drift, c1 + 10 * drift, 0.05, 0.5, height)
            for c2, c1, height in ((0.6, 6.0, 1.0), (1.4, 14.0, 0.8))
        )
        return Spectrum2D(ppm_f2, ppm_f1, image, Path("fine"))

    reference, shifted = fine_spectrum(0.0), fine_spectrum(0.02)
    common = dict(range_f2=(0.0, 2.0), range_f1=(0.0, 20.0), sigma_f2=0.03, sigma_f1=0.3)
    coarse = local_contrast_similarity(
        reference, shifted, step_f2=0.01, step_f1=0.1, **common
    )["similarity"]
    fine = local_contrast_similarity(
        reference, shifted, step_f2=0.005, step_f1=0.05, **common
    )["similarity"]
    assert abs(coarse - fine) < 1e-3


def test_local_contrast_handles_a_window_narrower_than_its_blur_kernel():
    spectrum = make_spectrum(REF, "spectrum")
    result = local_contrast_similarity(
        spectrum, spectrum, range_f2=(2.9, 3.1), range_f1=(39.0, 41.0), **KW
    )
    assert result["grid"] == [3, 5]
    assert 0.0 <= result["similarity"] <= 1.0


@pytest.mark.parametrize("name,value", [
    ("sigma_f2", 0.0), ("sigma_f1", -1.0),
    ("step_f2", 0.0), ("step_f1", -1.0),
    ("sigma_f2", np.nan), ("sigma_f1", np.inf),
    ("step_f2", np.nan), ("step_f1", np.inf),
])
def test_local_contrast_rejects_invalid_render_parameters(name, value):
    spectrum = make_spectrum(REF, "spectrum")
    with pytest.raises(ValueError, match=name):
        local_contrast_similarity(spectrum, spectrum, **{**KW, name: value})


def test_local_contrast_rejects_nonfinite_data_and_zero_feature_norm():
    spectrum = make_spectrum(REF, "spectrum")
    bad = Spectrum2D(
        spectrum.ppm_f2, spectrum.ppm_f1, spectrum.intensity.copy(), Path("bad")
    )
    bad.intensity[0, 0] = np.nan
    with pytest.raises(ValueError, match="spectrum_y"):
        local_contrast_similarity(spectrum, bad, **KW)
    with pytest.raises(ValueError, match="range_f2"):
        local_contrast_similarity(spectrum, spectrum, range_f2=(0.0, np.inf), **KW)

    assert local_contrast_similarity(
        spectrum, spectrum, sigma_f2=1e-4, sigma_f1=1e-3,
        step_f2=KW["step_f2"], step_f1=KW["step_f1"],
    )["similarity"] == 0.0


def test_cli_keeps_lcc_default_and_accepts_local_contrast():
    parser = build_parser()
    assert parser.parse_args(["x", "y"]).method == "lcc"
    assert parser.parse_args(["x", "y", "--method", "local-contrast"]).method == "local-contrast"
