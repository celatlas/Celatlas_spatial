def resolve_tissue_canvas_shape(tissue_bbox_metadata=None, h5_shape=None, coord_shape=None):
    """Resolve the full tissue canvas from current coordinate sources.

    Legacy ST inputs can store the vertical canvas extent in ``tissue_bbox``
    while the barcode coordinates extend the horizontal canvas beyond the
    nominal square chip width. Preserve that height when the coordinates are
    clearly expanding the width; otherwise an old square bbox is treated as
    stale metadata and ignored.
    """
    messages = []

    if h5_shape is not None:
        height, width = int(h5_shape[0]), int(h5_shape[1])
        source = 'barcodeToPos.h5'
        if coord_shape is not None:
            coord_h, coord_w = int(coord_shape[0]), int(coord_shape[1])
            if coord_h > height or coord_w > width:
                messages.append(
                    "Expanding barcodeToPos canvas to include FilterBarcodes "
                    f"coordinates: h5=({height}, {width}), coords=({coord_h}, {coord_w})"
                )
                height = max(height, coord_h)
                width = max(width, coord_w)
                source = 'barcodeToPos.h5 + FilterBarcodes'
    elif coord_shape is not None:
        height, width = int(coord_shape[0]), int(coord_shape[1])
        source = 'FilterBarcodes'
    else:
        raise FileNotFoundError(
            "Unable to resolve tissue canvas: need FilterBarcodes coordinates, "
            "or barcodeToPos.h5"
        )

    if tissue_bbox_metadata is not None:
        bbox = tuple(int(value) for value in tissue_bbox_metadata[:4])
        bbox_height, bbox_width = bbox[2], bbox[3]
        coordinate_width_expands_bbox = bbox_width < width
        if bbox_height > height and coordinate_width_expands_bbox:
            height = bbox_height
            source = f"tissue_bbox height + {source} width"
            messages.append(
                "Preserving tissue_bbox height while expanding canvas width from "
                "barcode coordinates"
            )
        else:
            messages.append("Ignoring tissue_bbox.csv for canvas resolution; using coordinate-derived canvas")

    return int(height), int(width), source, messages
