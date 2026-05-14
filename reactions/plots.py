"""Plotting utilities for the reactions library."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cycler

# Named color palette — consistent across all figures
COLORS = {
    "primary": "#1c4f8a",  # deep blue — main curve
    "secondary": "#e07b39",  # warm orange — contrast series
    "tertiary": "#2e8b57",  # sea green — third series
    "quaternary": "#8b5a9e",  # medium purple — fourth series
    "reference": "#888888",  # gray — dashed/dotted reference lines
}

# Figure sizes in inches
_SIZES = {
    "single": (6.0, 4.0),
    "wide": (9.0, 4.0),  # 1×2 subplots
    "panorama": (12.0, 3.5),  # 1×4 subplots
}


def setup_figure(
    nrows: int = 1,
    ncols: int = 1,
    **kwargs,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a styled figure. Drop-in replacement for plt.subplots().

    Parameters
    ----------
    nrows : int
        Number of subplot rows.
    ncols : int
        Number of subplot columns.
    **kwargs
        Forwarded to ``plt.subplots()``. ``figsize`` may be overridden.

    Returns
    -------
    fig, ax : Figure and Axes (or array of Axes for multi-panel figures).
    """
    if nrows == 1 and ncols == 1:
        key = "single"
    elif nrows == 1 and ncols == 2:
        key = "wide"
    elif nrows == 1 and ncols >= 4:
        key = "panorama"
    else:
        key = None

    figsize = kwargs.pop(
        "figsize",
        _SIZES[key] if key else (ncols * 6.0, nrows * 4.0),
    )
    return plt.subplots(nrows, ncols, figsize=figsize, **kwargs)


# Apply defaults globally on import
mpl.rcParams.update(
    {
        "lines.linewidth": 2.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": cycler(
            color=[
                COLORS["primary"],
                COLORS["secondary"],
                COLORS["tertiary"],
                COLORS["quaternary"],
            ]
        ),
        "legend.fontsize": 9,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        # 96 DPI = CSS web standard (1 CSS px = 1/96 in).
        # Wide figures (9 in) save at 864 px; typical column (~750 px) scales them
        # to ~87 %, so specified font sizes map predictably to rendered sizes.
        "savefig.dpi": 96,
    }
)
