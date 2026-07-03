from pathlib import Path

import numpy as np

from spectrum_similarity import Spectrum, read_bruker_1d, save_similarity_plot, spectrum_similarity


def gaussian(ppm, center, width, height=1.0):
    return height * np.exp(-0.5 * ((ppm - center) / width) ** 2)


def test_identical_spectra_have_similarity_one():
    ppm = np.linspace(0, 10, 4096)
    intensity = gaussian(ppm, 2.0, 0.03) + gaussian(ppm, 7.0, 0.05)
    spectrum = Spectrum(ppm=ppm, intensity=intensity, source=Path("synthetic"))

    result = spectrum_similarity(spectrum, spectrum, min_bin_width=0.4)

    assert result["similarity"] == 1.0


def test_shifted_related_spectrum_scores_higher_than_unrelated():
    ppm = np.linspace(0, 10, 4096)
    reference = Spectrum(
        ppm=ppm,
        intensity=gaussian(ppm, 2.0, 0.03) + gaussian(ppm, 7.0, 0.05),
        source=Path("reference"),
    )
    shifted = Spectrum(
        ppm=ppm,
        intensity=gaussian(ppm, 2.15, 0.03) + gaussian(ppm, 7.10, 0.05),
        source=Path("shifted"),
    )
    unrelated = Spectrum(
        ppm=ppm,
        intensity=gaussian(ppm, 4.0, 0.03) + gaussian(ppm, 9.0, 0.05),
        source=Path("unrelated"),
    )

    shifted_score = spectrum_similarity(reference, shifted, min_bin_width=0.4)["similarity"]
    unrelated_score = spectrum_similarity(reference, unrelated, min_bin_width=0.4)["similarity"]

    assert shifted_score > unrelated_score
    assert 0.0 <= unrelated_score <= 1.0
    assert 0.0 <= shifted_score <= 1.0


def test_reads_minimal_bruker_processed_1d(tmp_path):
    pdata = tmp_path / "sample" / "pdata" / "1"
    pdata.mkdir(parents=True)
    values = np.arange(8, dtype="<i4")
    values.tofile(pdata / "1r")
    (pdata / "procs").write_text(
        "\n".join(
            [
                "##$SI= 8",
                "##$OFFSET= 10",
                "##$SW_p= 5000",
                "##$SF= 500",
                "##$BYTORDP= 0",
                "##$DTYPP= 0",
                "##$NC_proc= 0",
            ]
        ),
        encoding="latin-1",
    )

    spectrum = read_bruker_1d(tmp_path / "sample")

    assert spectrum.source == pdata / "1r"
    assert spectrum.intensity.tolist() == list(range(8))
    assert spectrum.ppm[0] == 10.0
    assert spectrum.ppm[-1] == 1.25


def test_saves_similarity_plot(tmp_path):
    ppm = np.linspace(0, 10, 1024)
    reference = Spectrum(ppm=ppm, intensity=gaussian(ppm, 2.0, 0.03), source=Path("reference"))
    shifted = Spectrum(ppm=ppm, intensity=gaussian(ppm, 2.15, 0.03), source=Path("shifted"))
    result = spectrum_similarity(reference, shifted, min_bin_width=0.4)
    plot_path = tmp_path / "plot.png"

    saved_path = save_similarity_plot(result, reference, shifted, output=plot_path)

    assert saved_path == plot_path.resolve()
    assert plot_path.is_file()
    assert plot_path.stat().st_size > 0
