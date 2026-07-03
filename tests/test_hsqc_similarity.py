from pathlib import Path

import numpy as np

from hsqc_similarity import (
    Spectrum2D,
    hsqc_similarity,
    read_bruker_2d,
    save_similarity_plot_2d,
)


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


def test_identical_2d_spectra_have_similarity_one():
    spectrum = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0), (7.0, 120.0, 0.1, 1.5, 0.8)], "synthetic")

    result = hsqc_similarity(spectrum, spectrum, min_bin_width_f2=0.4, min_bin_width_f1=4.0)

    assert result["similarity"] == 1.0


def test_rotation_preserves_identity_and_ordering():
    reference = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0), (7.0, 120.0, 0.1, 1.5, 0.8)], "ref")
    shifted = make_spectrum([(3.1, 41.0, 0.1, 1.5, 1.0), (7.1, 121.0, 0.1, 1.5, 0.8)], "shifted")
    unrelated = make_spectrum([(5.0, 80.0, 0.1, 1.5, 1.0), (1.0, 20.0, 0.1, 1.5, 0.8)], "unrelated")

    assert hsqc_similarity(reference, reference, rotate_deg=45)["similarity"] == 1.0
    shifted_score = hsqc_similarity(reference, shifted, rotate_deg=45)["similarity"]
    unrelated_score = hsqc_similarity(reference, unrelated, rotate_deg=45)["similarity"]
    assert shifted_score > unrelated_score


def test_smoothing_preserves_identity():
    spectrum = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0), (7.0, 120.0, 0.1, 1.5, 0.8)], "synthetic")
    result = hsqc_similarity(spectrum, spectrum, smooth_f2=0.05, smooth_f1=0.5)
    assert result["similarity"] == 1.0


def test_related_2d_spectrum_scores_higher_than_unrelated():
    reference = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0), (7.0, 120.0, 0.1, 1.5, 0.8)], "ref")
    shifted = make_spectrum([(3.1, 41.0, 0.1, 1.5, 1.0), (7.1, 121.0, 0.1, 1.5, 0.8)], "shifted")
    unrelated = make_spectrum([(5.0, 80.0, 0.1, 1.5, 1.0), (1.0, 20.0, 0.1, 1.5, 0.8)], "unrelated")

    shifted_score = hsqc_similarity(reference, shifted)["similarity"]
    unrelated_score = hsqc_similarity(reference, unrelated)["similarity"]

    assert shifted_score > unrelated_score
    assert 0.0 <= unrelated_score <= 1.0
    assert 0.0 <= shifted_score <= 1.0


def _write_procs(path, si, offset, sw, sf, xdim):
    path.write_text(
        "\n".join(
            [
                f"##$SI= {si}",
                f"##$OFFSET= {offset}",
                f"##$SW_p= {sw}",
                f"##$SF= {sf}",
                f"##$XDIM= {xdim}",
                "##$BYTORDP= 0",
                "##$DTYPP= 0",
                "##$NC_proc= 0",
            ]
        ),
        encoding="latin-1",
    )


def test_reads_minimal_bruker_processed_2d(tmp_path):
    pdata = tmp_path / "sample" / "pdata" / "1"
    pdata.mkdir(parents=True)

    si_f2, si_f1 = 4, 4
    xdim_f2, xdim_f1 = 2, 2
    # Build a 4x4 matrix with unique values, store it in 2x2 submatrix (tile) order.
    matrix = np.arange(si_f1 * si_f2, dtype="<i4").reshape(si_f1, si_f2)
    tiles = []
    for tf1 in range(si_f1 // xdim_f1):
        for tf2 in range(si_f2 // xdim_f2):
            tile = matrix[
                tf1 * xdim_f1 : (tf1 + 1) * xdim_f1, tf2 * xdim_f2 : (tf2 + 1) * xdim_f2
            ]
            tiles.append(tile.ravel())
    np.concatenate(tiles).astype("<i4").tofile(pdata / "2rr")

    _write_procs(pdata / "procs", si_f2, offset=10, sw=5000, sf=500, xdim=xdim_f2)
    _write_procs(pdata / "proc2s", si_f1, offset=160, sw=20000, sf=125, xdim=xdim_f1)

    spectrum = read_bruker_2d(tmp_path / "sample")

    assert spectrum.source == pdata / "2rr"
    assert spectrum.intensity.shape == (si_f1, si_f2)
    # Submatrix reassembly must recover the original row-major matrix.
    np.testing.assert_array_equal(spectrum.intensity, matrix)
    assert spectrum.ppm_f2[0] == 10.0
    assert spectrum.ppm_f1[0] == 160.0


def test_saves_similarity_plot_2d(tmp_path):
    reference = make_spectrum([(3.0, 40.0, 0.1, 1.5, 1.0)], "ref")
    shifted = make_spectrum([(3.1, 41.0, 0.1, 1.5, 1.0)], "shifted")
    result = hsqc_similarity(reference, shifted)
    plot_path = tmp_path / "plot.png"

    saved_path = save_similarity_plot_2d(result, reference, shifted, output=plot_path)

    assert saved_path == plot_path.resolve()
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
