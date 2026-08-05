"""
Shared infrastructure for QC report generators
==============================================

This module is the **single source of truth** for the HTML shell, the panel
wrapper, the colour palettes, and the provenance carrier used by every
generator in ``trnaseq.qc``.

It exists so that the per-generator boilerplate stops growing as we add
panel-level features (narrative blocks, auto-findings, CSV exports, URL
anchors). Every generator in this package builds its body out of
:func:`render_panel` calls and finalises with :func:`render_html_shell`.

Public surface
--------------
* :class:`ReportContext`           — provenance & metadata for the footer.
* :func:`render_html_shell`        — emits ``<!DOCTYPE html>`` + ``<head>`` + ``<body>``.
* :func:`render_panel`             — wraps a single panel (h2, card, optional extras).
* :func:`export_panel_data`        — writes a per-panel CSV and returns a relative path.
* :func:`fig_to_div`               — uniform Plotly-figure-to-div conversion.
* :func:`resolve_trna_col`         — central column-name lookup (handles drift).
* :func:`host_filter`              — drop synthetic spike-ins; lifted from charge_report.
* ``STYLE_CSS``                    — shared stylesheet string.
* ``OKABE_ITO`` / ``GLASBEY_HEX``  — colour-blind-safe categorical palettes.

Notes
-----
T0 of the report-improvements PRP creates this module as a pure refactor:
``fig_div`` (the plot content) is **byte-identical** before/after — only the
surrounding HTML is centralised here. Later tasks (T1–T12) extend the same
helpers without touching the per-generator panel methods again.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import shutil
import socket
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import plotly.offline


# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------

#: 8-colour Okabe-Ito palette — colour-blind safe (Wong, Nature Methods 2011).
OKABE_ITO: list[str] = [
    '#E69F00',  # orange
    '#56B4E9',  # sky blue
    '#009E73',  # bluish green
    '#F0E442',  # yellow
    '#0072B2',  # blue
    '#D55E00',  # vermillion
    '#CC79A7',  # reddish purple
    '#000000',  # black
]

#: 24-colour palette derived from the Glasbey perceptually-distinct set
#: (pre-computed; avoids a colorcet dependency).
GLASBEY_HEX: list[str] = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78',
    '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7',
    '#dbdb8d', '#9edae5', '#393b79', '#637939', '#8c6d31', '#843c39',
]


def categorical_palette(n: int) -> list[str]:
    """Return ``n`` distinct hex colours.

    ≤8 categories → Okabe-Ito (best perceptual separation for small N).
    9–24 categories → GLASBEY_HEX.
    >24 categories → Glasbey + evenly-spaced HSV extras.
    """
    if n <= len(OKABE_ITO):
        return OKABE_ITO[:n]
    if n <= len(GLASBEY_HEX):
        return GLASBEY_HEX[:n]
    import colorsys
    extra = []
    deficit = n - len(GLASBEY_HEX)
    for i in range(deficit):
        h = ((i + 1) / (deficit + 1))
        r, g, b = colorsys.hls_to_rgb(h, 0.5, 0.7)
        extra.append(f'#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}')
    return GLASBEY_HEX + extra


# ---------------------------------------------------------------------------
# ReportContext — provenance carrier
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReportContext:
    """Provenance metadata threaded into every report footer.

    All fields are optional so the dataclass remains constructable in
    environments without git (tarball installs, unit tests).
    """
    project_dir: Optional[Path] = None
    pipeline_version: Optional[str] = None
    git_sha: Optional[str] = None
    git_branch: Optional[str] = None
    git_remote_url: Optional[str] = None
    config_path: Optional[Path] = None
    config_sha256: Optional[str] = None
    sample_sheet_path: Optional[Path] = None
    sample_sheet_sha256: Optional[str] = None
    reference_db: Optional[str] = None
    generated_at: _dt.datetime = field(default_factory=lambda: _dt.datetime.now(_dt.timezone.utc))
    runtime_seconds: Optional[float] = None
    host: str = field(default_factory=socket.gethostname)
    command: Optional[str] = field(default_factory=lambda: ' '.join(sys.argv) if sys.argv else None)

    # -- factory ----------------------------------------------------------

    @classmethod
    def discover(
        cls,
        project_dir: Optional[Path | str] = None,
        config_path: Optional[Path | str] = None,
        sample_sheet_path: Optional[Path | str] = None,
        reference_db: Optional[str] = None,
        pipeline_version: Optional[str] = None,
        runtime_seconds: Optional[float] = None,
        repo_root: Optional[Path | str] = None,
    ) -> 'ReportContext':
        """Build a context by sniffing git and hashing files.

        Tolerates every kind of failure (no git, missing files, permission
        errors) and falls back to ``None`` for any field it can't resolve.
        Never raises.
        """
        repo_root = Path(repo_root) if repo_root else _default_repo_root()
        return cls(
            project_dir=Path(project_dir) if project_dir else None,
            pipeline_version=pipeline_version or _read_pipeline_version(repo_root),
            git_sha=_git(['rev-parse', 'HEAD'], cwd=repo_root),
            git_branch=_git(['rev-parse', '--abbrev-ref', 'HEAD'], cwd=repo_root),
            git_remote_url=_git(['config', '--get', 'remote.origin.url'], cwd=repo_root),
            config_path=Path(config_path) if config_path else None,
            config_sha256=_sha256(config_path) if config_path else None,
            sample_sheet_path=Path(sample_sheet_path) if sample_sheet_path else None,
            sample_sheet_sha256=_sha256(sample_sheet_path) if sample_sheet_path else None,
            reference_db=reference_db,
            runtime_seconds=runtime_seconds,
        )

    # -- machine-readable export -----------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """JSON-safe dict (paths → strings, datetimes → ISO 8601)."""
        d = asdict(self)
        for k, v in list(d.items()):
            if isinstance(v, Path):
                d[k] = str(v)
            elif isinstance(v, _dt.datetime):
                d[k] = v.isoformat()
        return d

    def write_provenance_json(self, output_dir: Path | str) -> Path:
        """Write provenance.json next to the reports. Used by T9."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        path = output_dir / 'provenance.json'
        path.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True))
        return path

    # -- footer HTML ------------------------------------------------------

    def render_footer_html(self) -> str:
        """Return the HTML block that lives inside ``<div class="footer">``."""
        parts: list[str] = []

        # Pipeline line
        if self.git_sha:
            short = self.git_sha[:8]
            branch = self.git_branch or '?'
            commit_link = _github_commit_link(self.git_remote_url, self.git_sha)
            ver = f'tRNA-charge-seq {self.pipeline_version}' if self.pipeline_version else 'tRNA-charge-seq'
            if commit_link:
                parts.append(f'{ver} @ <a href="{commit_link}">{branch} ({short})</a>')
            else:
                parts.append(f'{ver} @ {branch} ({short})')
        else:
            parts.append('tRNA-charge-seq (git: unknown)')

        # Config + sample sheet
        if self.config_path:
            cfg_hash = (self.config_sha256 or '')[:12]
            parts.append(
                f'Config: <code title="sha256:{self.config_sha256}">'
                f'{self.config_path.name}</code> ({cfg_hash})'
            )
        if self.sample_sheet_path:
            ss_hash = (self.sample_sheet_sha256 or '')[:12]
            parts.append(
                f'Samples: <code title="sha256:{self.sample_sheet_sha256}">'
                f'{self.sample_sheet_path.name}</code> ({ss_hash})'
            )

        if self.reference_db:
            parts.append(f'Reference: <code>{self.reference_db}</code>')

        # Timestamp + host + runtime
        ts = self.generated_at.strftime('%Y-%m-%d %H:%M:%S %Z') or self.generated_at.isoformat()
        time_line = f'Generated: {ts} on <code>{self.host}</code>'
        if self.runtime_seconds is not None:
            time_line += f' in {self.runtime_seconds:.0f}s'
        parts.append(time_line)

        if self.command:
            parts.append(f'<code class="cmdline">{_html_escape(self.command)}</code>')

        return ' &nbsp;|&nbsp; '.join(parts[:-2]) + ('<br>' if len(parts) > 2 else '') + \
               ' &nbsp;|&nbsp; '.join(parts[-2:])


# ---------------------------------------------------------------------------
# Shared stylesheet
# ---------------------------------------------------------------------------

#: Single source of truth for the report stylesheet.
#:
#: Existing scrapers depend on the original class names (``.card``, ``.pass``,
#: ``.warn``, ``.fail``, ``.up``, ``.down``, ``.footer``) — they MUST keep
#: working. New classes (``.finding``, ``.dl-btn``, ``.panel-desc``,
#: ``.health-pill``) are additive for later tasks (T2/T9/T10).
STYLE_CSS: str = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       margin: 0; padding: 20px 40px; background: #f5f6fa; color: #2d3436; }
h1 { border-bottom: 3px solid var(--accent, #0984e3); padding-bottom: 10px; }
h2 { color: #2d3436; margin-top: 40px; scroll-margin-top: 20px; }
h2:target { background: #fff7d6; padding-left: 6px; border-left: 4px solid var(--accent, #0984e3); }
.card { background: white; border-radius: 8px; padding: 20px; margin: 20px 0;
        box-shadow: 0 2px 8px rgba(0,0,0,0.08); overflow-x: auto; }
.card:target { outline: 2px solid var(--accent, #0984e3); }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th { background: var(--accent, #0984e3); color: white; padding: 8px 12px; text-align: left;
     position: sticky; top: 0; }
td { padding: 6px 12px; border-bottom: 1px solid #dfe6e9; }
tr:hover { background: #f0f3f8; }
.pass { background-color: #00b89433; }
.warn { background-color: #fdcb6e55; }
.fail { background-color: #d6336c33; }
.up   { color: #d63031; font-weight: bold; }
.down { color: #0984e3; font-weight: bold; }
.footer { font-size: 11px; color: #636e72; margin-top: 40px; text-align: center;
          padding: 16px; border-top: 1px solid #dfe6e9; }
.footer code { background: #ecf0f1; padding: 1px 6px; border-radius: 3px;
               font-size: 11px; color: #2d3436; }
.footer code.cmdline { display: inline-block; max-width: 90%; overflow-x: auto;
                       white-space: nowrap; vertical-align: middle; }
.footer a { color: #0984e3; text-decoration: none; }
.footer a:hover { text-decoration: underline; }

/* T2 hooks — narrative + auto-findings */
.panel-desc { margin: 0 0 12px 0; font-size: 13px; color: #555; }
.panel-desc summary { cursor: pointer; color: #0984e3; font-weight: 500; }
.panel-desc p { margin: 6px 0 0 0; }
.finding { margin: 12px 0 0 0; padding: 8px 12px; border-radius: 4px;
           font-size: 13px; border-left: 4px solid; }
.finding.pass { background: #00b89422; border-color: #00b894; }
.finding.pass::before { content: "\\2713  "; }   /* check */
.finding.warn { background: #fdcb6e33; border-color: #f39c12; }
.finding.warn::before { content: "\\26A0  "; }   /* warning */
.finding.fail { background: #d6336c22; border-color: #d63031; }
.finding.fail::before { content: "\\2717  "; }   /* cross */

/* T10 hooks — CSV download */
.dl-btn { display: inline-block; margin: 12px 0 0 0; padding: 6px 14px;
          background: #ecf0f1; color: #2d3436; border-radius: 4px;
          font-size: 12px; text-decoration: none; border: 1px solid #cfd5d8; }
.dl-btn:hover { background: #dfe6e9; }
.dl-btn::before { content: "\\2B07  "; }   /* down arrow */

/* T1 hooks — index page tiles & health pills */
.health-pill { display: inline-block; padding: 2px 10px; border-radius: 12px;
               font-size: 11px; font-weight: 600; text-transform: uppercase;
               letter-spacing: 0.04em; }
.health-pill.pass { background: #00b894; color: white; }
.health-pill.warn { background: #f39c12; color: white; }
.health-pill.fail { background: #d63031; color: white; }
"""


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def render_html_shell(
    title: str,
    body: str,
    *,
    accent_color: str = '#0984e3',
    h1_text: Optional[str] = None,
    intro_html: str = '',
    context: Optional[ReportContext] = None,
    offline: bool = True,
    extra_head: str = '',
) -> str:
    """Return the full ``<!DOCTYPE html>`` … ``</html>`` page.

    Parameters
    ----------
    title : str
        Value for ``<title>``.
    body : str
        Concatenated panel HTML (typically ``'\\n'.join(panels)``).
    accent_color : str
        CSS colour applied to the H1 border + table header. Each report keeps
        its existing accent (charge #00b894, abundance/modification #6c5ce7,
        QC/fragment #0984e3).
    h1_text : str, optional
        Text inside ``<h1>``. Defaults to ``title``.
    intro_html : str
        Optional HTML inserted between ``<h1>`` and the panels (e.g. the
        ``<p>Project: ...</p>`` line in the QC report).
    context : ReportContext, optional
        Provenance info for the footer. If omitted, a minimal placeholder
        footer is emitted (matches pre-T0 behaviour).
    offline : bool
        If True (default), Plotly JS is inlined for offline use. T12 will
        flip the default to False (CDN with inline fallback).
    extra_head : str
        Extra HTML injected into ``<head>`` (used for Grid.js, Mermaid, etc.).
    """
    h1_text = h1_text if h1_text is not None else title

    if offline:
        plotly_script = f'<script>{plotly.offline.get_plotlyjs()}</script>'
    else:
        # T12 will fill in a CDN script + SRI hash + inline fallback.
        plotly_script = (
            '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" '
            'crossorigin="anonymous"></script>'
        )

    footer_html = context.render_footer_html() if context else 'Generated by trnaseq QC module'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{_html_escape(title)}</title>
{plotly_script}
<style>:root {{ --accent: {accent_color}; }}{STYLE_CSS}</style>
{extra_head}
</head>
<body>
<h1>{h1_text}</h1>
{intro_html}
{body}
<div class="footer">{footer_html}</div>
</body>
</html>"""


def render_panel(
    title: str,
    fig_div: str,
    *,
    anchor: Optional[str] = None,
    description: Optional[str] = None,
    finding: Optional[str] = None,
    severity: Optional[str] = None,
    csv_path: Optional[str] = None,
    extra_html: str = '',
) -> str:
    """Wrap a single panel.

    Output shape (slots empty unless arg provided)::

        <h2>{title}</h2>
        <div id="panel-{anchor}" class="card">
          <details class="panel-desc">…description…</details>   ← T2
          {fig_div}                                              ← unchanged
          <p class="finding {severity}">…</p>                    ← T2
          {extra_html}
          <a class="dl-btn" href="{csv_path}" download>…</a>     ← T10
        </div>

    For T0 every generator calls this with just ``title``, ``fig_div``, and
    ``anchor`` so the rendered output stays visually identical to the
    pre-refactor card (only ``id=`` is added to the card div).
    """
    id_attr = f' id="panel-{anchor}"' if anchor else ''

    pre_parts: list[str] = []
    if description:
        pre_parts.append(
            '<details class="panel-desc">'
            '<summary>How to read this</summary>'
            f'<p>{description}</p>'
            '</details>'
        )

    post_parts: list[str] = []
    if finding:
        sev = severity if severity in ('pass', 'warn', 'fail') else 'warn'
        post_parts.append(f'<p class="finding {sev}">{finding}</p>')
    if extra_html:
        post_parts.append(extra_html)
    if csv_path:
        post_parts.append(f'<a class="dl-btn" href="{csv_path}" download>Download CSV</a>')

    return (
        f'<h2>{title}</h2>'
        f'<div{id_attr} class="card">'
        f'{"".join(pre_parts)}{fig_div}{"".join(post_parts)}'
        f'</div>'
    )


def render_skipped_panel(title: str, reason: str, *, anchor: Optional[str] = None) -> str:
    """Render a placeholder card explaining why a panel was skipped.

    Per PRP §8 gotcha: "The summary code paths conditionally skip when a
    column is missing — but the report doesn't tell the user." This makes
    skips visible.
    """
    id_attr = f' id="panel-{anchor}"' if anchor else ''
    return (
        f'<h2>{title}</h2>'
        f'<div{id_attr} class="card" style="opacity:0.6">'
        f'<p style="margin:0;color:#636e72"><em>Panel skipped:</em> {reason}</p>'
        f'</div>'
    )


def fig_to_div(fig) -> str:
    """Convert a Plotly figure to an HTML div string (no plotly.js bundled).

    Centralised so all generators agree on how a figure is embedded.
    Equivalent to the per-generator ``_fig_to_div`` instance method.
    """
    return plotly.offline.plot(fig, output_type='div', include_plotlyjs=False)


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def export_panel_data(df, output_dir: Path | str, name: str) -> str:
    """Write ``<output_dir>/data/<name>.csv`` and return a relative href.

    The returned string is suitable for use as the ``csv_path=`` arg to
    :func:`render_panel`. The relative form keeps the report portable.
    """
    output_dir = Path(output_dir)
    data_dir = output_dir / 'data'
    data_dir.mkdir(parents=True, exist_ok=True)

    safe_name = name.replace('/', '_').replace('\\', '_')
    if not safe_name.endswith('.csv'):
        safe_name += '.csv'
    out_path = data_dir / safe_name
    df.to_csv(out_path, index=False)

    return f'data/{safe_name}'


# ---------------------------------------------------------------------------
# DataFrame helpers — single source of truth for column-name drift
# ---------------------------------------------------------------------------

#: Column-name candidates for the tRNA identifier (drift across generators).
_TRNA_COL_CANDIDATES: tuple[str, ...] = (
    'tRNA_annotation', 'tRNA_name', 'trna_name', 'feature', 'name',
)


def resolve_trna_col(df) -> Optional[str]:
    """Return the first present tRNA-identifier column name, or ``None``.

    Per PRP §8 gotcha: ``tRNA_annotation`` vs ``tRNA_name`` is handled
    inconsistently across generators — use this everywhere.
    """
    for candidate in _TRNA_COL_CANDIDATES:
        if candidate in df.columns:
            return candidate
    return None


def host_filter(df):
    """Filter ``df`` to host tRNAs (excluding synthetic spike-ins).

    Lifted from ``ChargeReportGenerator._host_filter`` per PRP §8 gotcha —
    abundance/fragment can now call this consistently.

    Returns the input unchanged when no source column is present.
    """
    if 'tRNA_source' in df.columns:
        return df[df['tRNA_source'] == 'host']
    if 'Syn_ctr' in df.columns:
        return df[~df['Syn_ctr']]
    return df


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def _git(args: list[str], cwd: Optional[Path]) -> Optional[str]:
    """Run ``git <args>`` quietly; return stripped stdout or ``None``."""
    if cwd is None or shutil.which('git') is None:
        return None
    try:
        res = subprocess.run(
            ['git', '-C', str(cwd)] + args,
            check=False, capture_output=True, text=True, timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if res.returncode != 0:
        return None
    out = res.stdout.strip()
    return out or None


def _sha256(path: Path | str | None) -> Optional[str]:
    """Stream-hash a file. Returns ``None`` on any failure."""
    if path is None:
        return None
    p = Path(path)
    if not p.is_file():
        return None
    try:
        h = hashlib.sha256()
        with p.open('rb') as fh:
            for chunk in iter(lambda: fh.read(65536), b''):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def _default_repo_root() -> Optional[Path]:
    """Best-effort repo root: walk up from this file until we hit a .git/."""
    here = Path(__file__).resolve()
    for parent in (here, *here.parents):
        if (parent / '.git').exists():
            return parent
    return None


def _read_pipeline_version(repo_root: Optional[Path]) -> Optional[str]:
    """Try to read a version from common locations; return ``None`` on miss."""
    if repo_root is None:
        return None
    candidates = [
        repo_root / 'trnaseq' / '__init__.py',
        repo_root / 'pyproject.toml',
        repo_root / 'VERSION',
    ]
    for c in candidates:
        if not c.is_file():
            continue
        try:
            text = c.read_text(encoding='utf-8', errors='ignore')
        except OSError:
            continue
        for line in text.splitlines():
            line = line.strip()
            if line.startswith('__version__') or line.startswith('version'):
                # __version__ = "1.2.3" or version = "1.2.3"
                for sep in ('=', ':'):
                    if sep in line:
                        rhs = line.split(sep, 1)[1].strip()
                        if rhs.startswith(('"', "'")):
                            return rhs.strip('"\' ,')
    return None


def _github_commit_link(remote_url: Optional[str], sha: Optional[str]) -> Optional[str]:
    """Turn ``git@github.com:org/repo.git`` (or https) + SHA into a commit URL."""
    if not remote_url or not sha:
        return None
    url = remote_url.strip()
    if url.endswith('.git'):
        url = url[:-4]
    if url.startswith('git@github.com:'):
        url = 'https://github.com/' + url[len('git@github.com:'):]
    elif url.startswith('ssh://git@github.com/'):
        url = 'https://github.com/' + url[len('ssh://git@github.com/'):]
    if not url.startswith('http'):
        return None
    return f'{url}/commit/{sha}'


def _html_escape(s: Any) -> str:
    """Tiny HTML escape (avoid pulling in html stdlib at call sites)."""
    text = str(s)
    return (text.replace('&', '&amp;')
                .replace('<', '&lt;')
                .replace('>', '&gt;')
                .replace('"', '&quot;'))


__all__ = [
    'ReportContext',
    'STYLE_CSS',
    'OKABE_ITO',
    'GLASBEY_HEX',
    'categorical_palette',
    'render_html_shell',
    'render_panel',
    'render_skipped_panel',
    'fig_to_div',
    'export_panel_data',
    'resolve_trna_col',
    'host_filter',
]
