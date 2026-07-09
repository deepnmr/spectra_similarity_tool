import json
import tempfile
from pathlib import Path

import pytest

from bench_13c import _assert_distinct, load_peaklist
from hsqc_lcc import lcc_similarity


def _write_peaklist(peaks):
    """peaks: list of (1H_ppm, 13C_ppm, intensity) -> temp JSON path in the simpleNMR shape."""
    idx = [str(i + 1) for i in range(len(peaks))]
    hsqc = {
        "f2_ppm": {k: p[0] for k, p in zip(idx, peaks)},
        "f1_ppm": {k: p[1] for k, p in zip(idx, peaks)},
        "intensity": {k: p[2] for k, p in zip(idx, peaks)},
        "signaltype": {k: "Compound" for k in idx},
    }
    f = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump({"hsqc": hsqc}, f)
    f.close()
    return Path(f.name)


def test_loader_places_peaks_and_ignores_sign():
    # A negative intensity (edited-HSQC CH2) must still render a peak (uses |intensity|).
    path = _write_peaklist([(3.26, 71.1, 2.4), (1.86, 44.9, -1.9)])
    spec = load_peaklist(path)
    # both peaks present as positive intensity in the rendered image
    assert spec.intensity.max() > 0
    i_h = int(round(3.26 / 0.01))
    i_c = int(round(71.1 / 0.10))
    assert spec.intensity[i_c, i_h] > 0


def test_loader_self_similarity_is_one():
    path = _write_peaklist([(3.26, 71.1, 2.4), (0.9, 22.0, 1.0), (7.2, 128.0, 0.5)])
    spec = load_peaklist(path)
    assert lcc_similarity(spec, spec, sigma_f2=0.05, sigma_f1=0.5, step_f2=0.02, step_f1=0.2,
                          range_f2=(0.0, 10.0), range_f1=(0.0, 165.0))["similarity"] == 1.0


def test_duplicate_same_pair_is_rejected():
    # A same-compound pair must be two distinct recordings; byte-identical files (the
    # olivetol bug) render to an identical image and must be rejected, not scored as 1.00.
    spec = load_peaklist(_write_peaklist([(3.26, 71.1, 2.4), (0.9, 22.0, 1.0)]))
    with pytest.raises(ValueError):
        _assert_distinct({"a": spec, "b": spec}, [("a", "b")])
