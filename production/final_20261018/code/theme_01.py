"""Shared rendering style extracted from the user's 01 figure package.

This module only controls appearance. It never reads clinical inputs, changes
scores, recomputes tests, or changes gene-set membership. Font files are not bundled.
"""
from __future__ import annotations

import os
import json
import warnings
from pathlib import Path
from typing import Any

import matplotlib as mpl
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(os.environ['ICI_PRODUCTION_ROOT'])


def load_style(path: str | Path | None = None) -> dict[str, Any]:
    path = Path(path) if path else ROOT / 'config' / 'style_01.json'
    return json.loads(path.read_text(encoding='utf-8'))


def resolve_font(style: dict[str, Any]) -> tuple[str, bool]:
    requested = style['production_targets']['font']
    available = {f.name for f in font_manager.fontManager.ttflist}
    for family in [requested, *style['production_targets']['fallback_order']]:
        if family in available:
            return family, family != requested
    raise RuntimeError('No configured font is installed. Install an authorized font locally; do not extract/share font files.')


def apply_style(path: str | Path | None = None) -> dict[str, Any]:
    style = load_style(path)
    target = style['production_targets']
    palette = style['palette']
    font, fallback = resolve_font(style)
    if fallback:
        warnings.warn(f"Requested Arial is unavailable; using {font}. Record this in the final font report.", stacklevel=2)
    mpl.rcParams.update({
        'font.family': font,
        'font.size': target['body_pt'],
        'axes.labelsize': target['axis_pt'],
        'axes.titlesize': target['panel_title_pt'],
        'axes.titleweight': 'bold',
        'axes.titlelocation': 'left',
        'xtick.labelsize': target['tick_pt'],
        'ytick.labelsize': target['tick_pt'],
        'legend.fontsize': target['tick_pt'],
        'figure.facecolor': palette['white'],
        'axes.facecolor': palette['white'],
        'savefig.facecolor': palette['white'],
        'text.color': palette['ink'],
        'axes.labelcolor': palette['ink'],
        'xtick.color': palette['ink'],
        'ytick.color': palette['ink'],
        'axes.edgecolor': palette['muted'],
        'axes.linewidth': target['axes_lw'],
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': False,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': target['axes_lw'],
        'ytick.major.width': target['axes_lw'],
        'lines.linewidth': target['data_line_lw'],
        'lines.markersize': target['point_size_pt'],
        'legend.frameon': False,
        'pdf.fonttype': target['pdf_fonttype'],
        'ps.fonttype': 42,
        'svg.fonttype': target['svg_fonttype'],
        'savefig.transparent': False,
    })
    return {'font_requested': target['font'], 'font_used': font, 'fallback': fallback, 'template': style['template']}


def signed_cmap(path: str | Path | None = None) -> LinearSegmentedColormap:
    """Colour represents the signed plotted variable, not clinical benefit."""
    hm = load_style(path)['heatmap']
    cmap = LinearSegmentedColormap.from_list('ici_01_signed', [hm['negative'], hm['zero'], hm['positive']])
    cmap.set_bad(hm['ne'])
    return cmap


def mm_to_inches(mm: float) -> float:
    if mm <= 0:
        raise ValueError('A physical dimension must be positive.')
    return mm / 25.4


def finish_axis(ax, *, zero: str | None = None) -> None:
    """Apply shared axes and optional labelled-data zero reference."""
    st = load_style()
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=2.5, width=0.5)
    if zero in ('x', 'both'):
        ax.axvline(0, color=st['palette']['muted'], linewidth=0.5, linestyle=(0, (3, 3)), zorder=0)
    if zero in ('y', 'both'):
        ax.axhline(0, color=st['palette']['muted'], linewidth=0.5, linestyle=(0, (3, 3)), zorder=0)


def save_vector(fig, stem: str | Path, *, png: bool = True) -> list[str]:
    """Export at the figure's physical size; do not silently tight-crop dimensions."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    exports = []
    for ext in ('pdf', 'svg'):
        path = stem.with_suffix('.' + ext)
        fig.savefig(path, bbox_inches=None)
        exports.append(str(path))
    if png:
        path = stem.with_suffix('.png')
        fig.savefig(path, dpi=180, bbox_inches=None)
        exports.append(str(path))
    return exports
