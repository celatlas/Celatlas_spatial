"""Method/module registry for Celatlas runner planning."""

from __future__ import annotations

from dataclasses import dataclass


class MethodError(ValueError):
    """Raised when a method cannot be resolved."""


@dataclass(frozen=True)
class MethodModule:
    """Workflow-independent behavior for a sample method.

    ``backend_method`` is the value passed to binSegment and report steps.
    Public methods are kept canonical; compatibility aliases are normalized
    before this registry is queried.
    """

    name: str
    backend_method: str
    uses_gene_expression: bool = False
    requires_tissue_tif: bool = False
    requires_he_image: bool = False

    def extend_binsegment_command(self, command: list[str], *, sample: str, input_dir, he_image=None, values=None) -> None:
        values = values or {}
        command.extend(["--method", self.backend_method])

        if self.requires_tissue_tif:
            command.extend(
                [
                    "--tif",
                    str(input_dir / f"{sample}.tif"),
                    "--ssdna-threshold-scale",
                    str(values["ssdna_threshold_scale"]),
                    "--ssdna-mask-expand-pixels",
                    str(values["ssdna_mask_expand_pixels"]),
                    "--ssdna-min-hole-area",
                    str(values.get("ssdna_min_hole_area", 5000)),
                ]
            )

        if self.uses_gene_expression:
            command.extend(
                [
                    "--gem-bin-size",
                    str(values["gem_bin_size"]),
                    "--umi-min-threshold",
                    str(values["umi_min_threshold"]),
                    "--enhance-params",
                    str(values["enhance_params"]),
                ]
            )

        if self.requires_he_image:
            if he_image is None:
                raise MethodError(f"Method '{self.name}' requires an H&E image path.")
            command.extend(["--tif", str(he_image), "--registration-type", str(values["registration_type"])])


_MODULES = {
    "gene_expr": MethodModule("gene_expr", "gene_expr", uses_gene_expression=True),
    "ssDNA": MethodModule("ssDNA", "ssDNA", requires_tissue_tif=True),
    "HE": MethodModule("HE", "HE", uses_gene_expression=True, requires_he_image=True),
}

_ALIASES = {
    "gene_expr": "gene_expr",
    "gene-expression": "gene_expr",
    "gene_expression": "gene_expr",
    "expression": "gene_expr",
    "image": "ssDNA",
    "tif": "ssDNA",
    "tiff": "ssDNA",
    "ssdna": "ssDNA",
    "ss_dna": "ssDNA",
    "ss-dna": "ssDNA",
    "he": "HE",
    "h&e": "HE",
    "h_and_e": "HE",
    "h-e": "HE",
}


def get_method_module(method: str) -> MethodModule:
    key = normalize_method(method)
    try:
        return _MODULES[key]
    except KeyError as exc:
        expected = ", ".join(sorted(_MODULES))
        raise MethodError(f"Unknown method '{method}'. Expected one of: {expected}.") from exc


def normalize_method(method: str) -> str:
    key = (method or "").strip()
    alias_key = key.lower()
    return _ALIASES.get(alias_key, key)


def backend_method_for(method: str) -> str:
    return get_method_module(method).backend_method


def requires_tissue_tif(method: str) -> bool:
    return get_method_module(method).requires_tissue_tif


def requires_he_image(method: str) -> bool:
    return get_method_module(method).requires_he_image
