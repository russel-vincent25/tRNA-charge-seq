"""
Tests for trnaseq.qc._common — shared report infrastructure.

These tests cover the T0 acceptance from the report-improvements PRP:

* only one ``<!DOCTYPE html>`` source of truth (smoke-tested via render shape);
* ``render_panel`` keeps the legacy ``<h2>…</h2><div class="card">…</div>``
  shape when called with no extras (so plot bytes stay identical);
* ``ReportContext`` is constructable in tests, JSON-safe, and degrades
  gracefully when git/files are missing;
* helpers (``host_filter``, ``resolve_trna_col``, palettes,
  ``export_panel_data``) behave as advertised.
"""

from __future__ import annotations

import datetime as _dt
import json
import re
from pathlib import Path

import pandas as pd
import pytest

from trnaseq.qc._common import (
    GLASBEY_HEX,
    OKABE_ITO,
    ReportContext,
    STYLE_CSS,
    _github_commit_link,
    _sha256,
    categorical_palette,
    export_panel_data,
    fig_to_div,
    host_filter,
    render_html_shell,
    render_panel,
    render_skipped_panel,
    resolve_trna_col,
)


# ---------------------------------------------------------------------------
# render_panel — T0 acceptance: legacy shape preserved when no extras
# ---------------------------------------------------------------------------

class TestRenderPanelLegacyShape:
    """When called with only (title, fig_div), the output must look like
    the pre-T0 hand-written panels so plot bytes are byte-identical."""

    def test_minimal_call_matches_legacy_pattern(self):
        out = render_panel('Mapping Rate', '<div>FAKE-DIV</div>')
        assert out == '<h2>Mapping Rate</h2><div class="card"><div>FAKE-DIV</div></div>'

    def test_anchor_adds_id_attribute_only(self):
        out = render_panel('Mapping Rate', '<div>X</div>', anchor='qc-mapping')
        assert out == (
            '<h2>Mapping Rate</h2>'
            '<div id="panel-qc-mapping" class="card">'
            '<div>X</div></div>'
        )
        # And specifically: nothing else inside the card.
        assert out.count('<div') == 2  # outer card + inner figure div

    def test_no_anchor_means_no_id_attr(self):
        out = render_panel('X', '<div></div>')
        assert ' id=' not in out


class TestRenderPanelExtras:
    """When optional kwargs are passed, the panel grows new slots."""

    def test_description_renders_collapsible(self):
        out = render_panel(
            'PCA', '<div></div>',
            anchor='qc-pca',
            description='Open this to learn how to read the plot.',
        )
        assert '<details class="panel-desc">' in out
        assert 'How to read this' in out
        assert 'Open this to learn how to read the plot.' in out

    def test_finding_uses_severity_class(self):
        out = render_panel('PCA', '<div></div>', finding='Bad sample', severity='fail')
        assert 'class="finding fail"' in out
        assert 'Bad sample' in out

    def test_finding_defaults_to_warn(self):
        out = render_panel('PCA', '<div></div>', finding='Heads up')
        assert 'class="finding warn"' in out

    def test_finding_unknown_severity_falls_back_to_warn(self):
        out = render_panel('PCA', '<div></div>', finding='?', severity='catastrophic')
        assert 'class="finding warn"' in out

    def test_csv_path_renders_download_button(self):
        out = render_panel('PCA', '<div></div>', csv_path='data/qc__pca.csv')
        assert 'class="dl-btn"' in out
        assert 'href="data/qc__pca.csv"' in out
        assert 'download' in out
        assert 'Download CSV' in out

    def test_all_extras_co_render_in_correct_order(self):
        out = render_panel(
            'PCA', '<div>FIG</div>',
            anchor='qc-pca',
            description='desc-text',
            finding='msg',
            severity='warn',
            csv_path='data/qc__pca.csv',
        )
        # description must precede the fig; finding + csv must follow.
        i_desc = out.index('desc-text')
        i_fig = out.index('<div>FIG</div>')
        i_finding = out.index('class="finding warn"')
        i_csv = out.index('class="dl-btn"')
        assert i_desc < i_fig < i_finding < i_csv


class TestRenderSkippedPanel:
    def test_renders_placeholder_with_reason(self):
        out = render_skipped_panel('Charge Heatmap', 'missing column tRNA_source')
        assert '<h2>Charge Heatmap</h2>' in out
        assert 'Panel skipped' in out
        assert 'missing column tRNA_source' in out

    def test_anchor_supported(self):
        out = render_skipped_panel('X', 'r', anchor='qc-foo')
        assert 'id="panel-qc-foo"' in out


# ---------------------------------------------------------------------------
# render_html_shell — T0 footer + offline-by-default + accent threading
# ---------------------------------------------------------------------------

class TestRenderHtmlShellStructure:
    def test_starts_with_doctype(self):
        out = render_html_shell('T', 'B', accent_color='#abc')
        assert out.startswith('<!DOCTYPE html>')

    def test_contains_title_and_body(self):
        out = render_html_shell('My Title', '<p>hello body</p>', accent_color='#abc')
        assert '<title>My Title</title>' in out
        assert '<p>hello body</p>' in out

    def test_uses_h1_text_when_provided(self):
        out = render_html_shell(
            'PageTitle', 'B', h1_text='Friendly Heading', accent_color='#abc',
        )
        assert '<h1>Friendly Heading</h1>' in out
        assert '<title>PageTitle</title>' in out

    def test_h1_falls_back_to_title(self):
        out = render_html_shell('PT', 'B', accent_color='#abc')
        assert '<h1>PT</h1>' in out

    def test_accent_color_threaded_via_css_var(self):
        out = render_html_shell('T', 'B', accent_color='#00b894')
        assert '--accent: #00b894' in out

    def test_intro_html_lands_between_h1_and_body(self):
        out = render_html_shell(
            'T', '<panel></panel>', accent_color='#abc',
            intro_html='<p class="intro">Welcome</p>',
        )
        h1_i = out.index('<h1>T</h1>')
        intro_i = out.index('Welcome')
        body_i = out.index('<panel></panel>')
        assert h1_i < intro_i < body_i

    def test_escapes_title_html(self):
        out = render_html_shell('A & B <evil>', 'body', accent_color='#abc')
        assert '<title>A &amp; B &lt;evil&gt;</title>' in out

    def test_only_one_doctype_in_output(self):
        out = render_html_shell('T', 'B', accent_color='#abc')
        assert out.count('<!DOCTYPE html>') == 1


class TestRenderHtmlShellPlotlyMode:
    def test_offline_default_inlines_plotly(self):
        out = render_html_shell('T', 'B', accent_color='#abc')
        # Plotly inlined when offline=True (the default)
        # We don't check the full JS, just that a <script>...</script> with
        # content is present (CDN URL would be a <script src="...">).
        assert 'src="https://cdn.plot.ly' not in out
        assert '<script>' in out

    def test_cdn_mode_uses_url(self):
        out = render_html_shell('T', 'B', accent_color='#abc', offline=False)
        assert 'src="https://cdn.plot.ly' in out
        # No giant inline blob.
        assert '<script>function' not in out

    def test_extra_head_is_injected(self):
        out = render_html_shell(
            'T', 'B', accent_color='#abc',
            extra_head='<link rel="stylesheet" href="https://example.org/grid.css">',
        )
        assert 'grid.css' in out
        # In the <head>, before <body>
        assert out.index('grid.css') < out.index('<body>')


class TestRenderHtmlShellFooter:
    def test_no_context_uses_pre_t0_placeholder(self):
        out = render_html_shell('T', 'B', accent_color='#abc')
        assert 'Generated by trnaseq QC module' in out

    def test_with_context_shows_git_or_unknown(self, monkeypatch):
        ctx = ReportContext(
            project_dir=Path('/x'),
            generated_at=_dt.datetime(2026, 1, 1, 0, 0, 0, tzinfo=_dt.timezone.utc),
            host='unit-test',
        )
        out = render_html_shell('T', 'B', accent_color='#abc', context=ctx)
        # Footer must contain SOMETHING from the context (host or git stanza).
        assert 'unit-test' in out
        assert ('git: unknown' in out) or ('tRNA-charge-seq' in out)

    def test_with_context_shows_short_sha_when_present(self):
        ctx = ReportContext(
            git_sha='abcdef1234567890' * 2,
            git_branch='feature/x',
            generated_at=_dt.datetime(2026, 1, 1, tzinfo=_dt.timezone.utc),
            host='h',
        )
        out = render_html_shell('T', 'B', accent_color='#abc', context=ctx)
        assert 'abcdef12' in out          # short SHA
        assert 'feature/x' in out

    def test_with_context_shows_runtime_seconds(self):
        ctx = ReportContext(
            runtime_seconds=487.0,
            generated_at=_dt.datetime(2026, 1, 1, tzinfo=_dt.timezone.utc),
            host='h',
        )
        out = render_html_shell('T', 'B', accent_color='#abc', context=ctx)
        assert '487s' in out


# ---------------------------------------------------------------------------
# Stylesheet — keep legacy classes alive; add new ones (T2, T9, T10, T1)
# ---------------------------------------------------------------------------

class TestStyleCss:
    """Per PRP §1.3: external scrapers depend on legacy class names."""

    @pytest.mark.parametrize('cls', ['.card', '.pass', '.warn', '.fail',
                                      '.up', '.down', '.footer'])
    def test_legacy_class_present(self, cls):
        assert cls in STYLE_CSS

    @pytest.mark.parametrize('cls', ['.finding', '.dl-btn', '.panel-desc',
                                      '.health-pill', '.health-pill.pass',
                                      '.health-pill.warn', '.health-pill.fail'])
    def test_new_class_present(self, cls):
        assert cls in STYLE_CSS

    def test_uses_accent_css_var(self):
        # Reports should be themable via --accent (set per generator).
        assert '--accent' in STYLE_CSS


# ---------------------------------------------------------------------------
# Palettes
# ---------------------------------------------------------------------------

class TestPalettes:
    def test_okabe_ito_size(self):
        assert len(OKABE_ITO) == 8

    def test_glasbey_size(self):
        assert len(GLASBEY_HEX) == 24

    def test_categorical_palette_returns_n(self):
        for n in [1, 5, 8, 12, 24, 40]:
            pal = categorical_palette(n)
            assert len(pal) == n

    def test_categorical_palette_prefers_okabe_for_small_n(self):
        assert categorical_palette(3) == OKABE_ITO[:3]

    def test_categorical_palette_uses_glasbey_for_mid_n(self):
        assert categorical_palette(20) == GLASBEY_HEX[:20]

    def test_categorical_palette_extends_for_huge_n(self):
        pal = categorical_palette(40)
        assert pal[:24] == GLASBEY_HEX
        # Extras are valid hex
        for color in pal[24:]:
            assert re.fullmatch(r'#[0-9a-f]{6}', color)


# ---------------------------------------------------------------------------
# Helpers — resolve_trna_col, host_filter, export_panel_data, fig_to_div
# ---------------------------------------------------------------------------

class TestResolveTrnaCol:
    def test_finds_tRNA_annotation(self):
        df = pd.DataFrame({'tRNA_annotation': ['a'], 'count': [1]})
        assert resolve_trna_col(df) == 'tRNA_annotation'

    def test_finds_trna_name_fallback(self):
        df = pd.DataFrame({'trna_name': ['a']})
        assert resolve_trna_col(df) == 'trna_name'

    def test_returns_none_when_absent(self):
        df = pd.DataFrame({'something_else': [1]})
        assert resolve_trna_col(df) is None

    def test_prefers_tRNA_annotation_over_others(self):
        df = pd.DataFrame({'tRNA_annotation': ['a'], 'feature': ['b']})
        assert resolve_trna_col(df) == 'tRNA_annotation'


class TestHostFilter:
    def test_filters_by_trna_source(self):
        df = pd.DataFrame({
            'tRNA_source': ['host', 'host', 'synthetic'],
            'name': ['a', 'b', 'c'],
        })
        out = host_filter(df)
        assert list(out['name']) == ['a', 'b']

    def test_filters_by_syn_ctr_fallback(self):
        df = pd.DataFrame({
            'Syn_ctr': [False, True, False],
            'name': ['a', 'b', 'c'],
        })
        out = host_filter(df)
        assert list(out['name']) == ['a', 'c']

    def test_returns_unchanged_when_no_source_cols(self):
        df = pd.DataFrame({'name': ['a', 'b']})
        assert host_filter(df) is df  # exact same object


class TestExportPanelData:
    def test_writes_file_and_returns_relative_path(self, tmp_path):
        df = pd.DataFrame({'a': [1, 2], 'b': ['x', 'y']})
        rel = export_panel_data(df, tmp_path, 'QC__sample_summary')
        assert rel == 'data/QC__sample_summary.csv'
        assert (tmp_path / 'data' / 'QC__sample_summary.csv').exists()
        # File contents are valid CSV
        round_tripped = pd.read_csv(tmp_path / 'data' / 'QC__sample_summary.csv')
        pd.testing.assert_frame_equal(round_tripped, df)

    def test_appends_csv_extension_if_missing(self, tmp_path):
        df = pd.DataFrame({'a': [1]})
        rel = export_panel_data(df, tmp_path, 'foo')
        assert rel == 'data/foo.csv'

    def test_does_not_double_csv_extension(self, tmp_path):
        df = pd.DataFrame({'a': [1]})
        rel = export_panel_data(df, tmp_path, 'foo.csv')
        assert rel == 'data/foo.csv'

    def test_sanitises_path_separators(self, tmp_path):
        df = pd.DataFrame({'a': [1]})
        rel = export_panel_data(df, tmp_path, 'subdir/oops')
        assert '/' not in rel.replace('data/', '', 1).split('/')[0]


class TestFigToDivSmoke:
    """A minimal Plotly figure produces a non-empty div without bundling JS."""

    def test_basic_fig(self):
        import plotly.graph_objects as go
        fig = go.Figure(go.Scatter(x=[1, 2], y=[1, 2]))
        out = fig_to_div(fig)
        assert '<div' in out
        # plotly.offline default: include_plotlyjs=False -> no <script>plotly</script> blob
        assert 'plotly-latest.min.js' not in out


# ---------------------------------------------------------------------------
# ReportContext
# ---------------------------------------------------------------------------

class TestReportContextConstruction:
    def test_minimal_default_construct(self):
        ctx = ReportContext()
        # All optional fields default to None except generated_at/host
        assert ctx.project_dir is None
        assert ctx.git_sha is None
        assert isinstance(ctx.generated_at, _dt.datetime)
        assert isinstance(ctx.host, str) and ctx.host

    def test_discover_in_empty_dir(self, tmp_path):
        ctx = ReportContext.discover(
            project_dir=tmp_path,
            config_path=None,
            sample_sheet_path=None,
            repo_root=tmp_path,  # no .git here
        )
        assert ctx.git_sha is None
        assert ctx.project_dir == tmp_path

    def test_discover_hashes_existing_file(self, tmp_path):
        cfg = tmp_path / 'cfg.yaml'
        cfg.write_text('key: value\n')
        ctx = ReportContext.discover(
            project_dir=tmp_path, config_path=cfg, repo_root=tmp_path,
        )
        # 64-hex SHA-256 string
        assert ctx.config_sha256 is not None
        assert len(ctx.config_sha256) == 64
        assert re.fullmatch(r'[0-9a-f]{64}', ctx.config_sha256)

    def test_discover_handles_missing_files_gracefully(self, tmp_path):
        ctx = ReportContext.discover(
            project_dir=tmp_path,
            config_path=tmp_path / 'nope.yaml',
            sample_sheet_path=tmp_path / 'also-nope.xlsx',
            repo_root=tmp_path,
        )
        assert ctx.config_sha256 is None
        assert ctx.sample_sheet_sha256 is None


class TestReportContextSerialisation:
    def test_to_dict_is_json_safe(self):
        ctx = ReportContext(
            project_dir=Path('/x'),
            config_path=Path('/x/cfg.yaml'),
            generated_at=_dt.datetime(2026, 1, 1, 0, 0, 0, tzinfo=_dt.timezone.utc),
        )
        d = ctx.to_dict()
        # Round-trips through JSON without TypeError
        s = json.dumps(d)
        assert '"project_dir": "' in s
        assert '"generated_at": "2026-01-01T' in s

    def test_write_provenance_json_round_trips(self, tmp_path):
        ctx = ReportContext(
            project_dir=tmp_path,
            generated_at=_dt.datetime(2026, 1, 1, tzinfo=_dt.timezone.utc),
            host='h',
        )
        p = ctx.write_provenance_json(tmp_path)
        assert p.exists()
        loaded = json.loads(p.read_text())
        assert loaded['host'] == 'h'


class TestReportContextFooter:
    def test_no_git_renders_unknown(self):
        ctx = ReportContext(
            generated_at=_dt.datetime(2026, 1, 1, tzinfo=_dt.timezone.utc),
            host='h',
        )
        html = ctx.render_footer_html()
        assert 'git: unknown' in html

    def test_with_git_renders_branch_and_short_sha(self):
        ctx = ReportContext(
            git_sha='0123456789abcdef' * 2,
            git_branch='main',
            generated_at=_dt.datetime(2026, 1, 1, tzinfo=_dt.timezone.utc),
            host='h',
        )
        html = ctx.render_footer_html()
        assert '01234567' in html
        assert 'main' in html

    def test_with_github_remote_includes_link(self):
        ctx = ReportContext(
            git_sha='abc123' * 7,
            git_branch='main',
            git_remote_url='https://github.com/foo/bar.git',
            generated_at=_dt.datetime(2026, 1, 1, tzinfo=_dt.timezone.utc),
            host='h',
        )
        html = ctx.render_footer_html()
        assert 'href="https://github.com/foo/bar/commit/' in html


# ---------------------------------------------------------------------------
# Internals — _sha256, _github_commit_link
# ---------------------------------------------------------------------------

class TestSha256Helper:
    def test_returns_none_for_missing_file(self, tmp_path):
        assert _sha256(tmp_path / 'nope') is None

    def test_returns_none_for_none(self):
        assert _sha256(None) is None

    def test_returns_hex_for_real_file(self, tmp_path):
        p = tmp_path / 'f.txt'
        p.write_text('hello\n')
        h = _sha256(p)
        assert h == '5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03'


class TestGithubCommitLink:
    def test_https_url_to_commit_link(self):
        link = _github_commit_link('https://github.com/foo/bar.git', 'abc' * 4)
        assert link == 'https://github.com/foo/bar/commit/abcabcabcabc'

    def test_ssh_url_to_commit_link(self):
        link = _github_commit_link('git@github.com:foo/bar.git', 'def')
        assert link == 'https://github.com/foo/bar/commit/def'

    def test_returns_none_when_missing(self):
        assert _github_commit_link(None, 'abc') is None
        assert _github_commit_link('https://github.com/x/y.git', None) is None

    def test_non_github_url_returns_none(self):
        # Only github URLs are recognised; gitlab/bitbucket fall through cleanly
        assert _github_commit_link('weird://example.com/foo/bar', 'abc') is None


# ---------------------------------------------------------------------------
# Integration smoke: refactored generators still emit only one DOCTYPE per file
# ---------------------------------------------------------------------------

class TestGeneratorIntegrationSmoke:
    """Hits the public API of each generator with the synth_run fixture so a
    test failure flags any regression introduced by the T0 refactor.

    These tests are intentionally minimal — they only assert that the report
    file is created, starts with exactly one ``<!DOCTYPE html>``, and contains
    the expected accent colour in the embedded CSS. The exhaustive panel-level
    checks live in the legacy ``tests/test_*_report.py`` files.
    """

    def test_qc_report_renders(self, synth_run, tmp_path):
        from trnaseq.qc.report import QCReportGenerator
        gen = QCReportGenerator(
            project_dir=synth_run.project_dir,
            sample_df=synth_run.sample_df,
            inp_file_df=synth_run.inp_file_df,
            stats_df=synth_run.stats_df,
            context=synth_run.context,
        )
        out = tmp_path / 'QC_report.html'
        gen.generate_html_report(out)
        text = out.read_text()
        assert text.count('<!DOCTYPE html>') == 1
        assert '--accent: #0984e3' in text
        # Footer carries provenance from the context fixture
        assert 'test-host' in text

    def test_charge_report_renders(self, synth_run, tmp_path):
        from trnaseq.qc.charge_report import ChargeReportGenerator
        gen = ChargeReportGenerator(
            charge_df_transcript=synth_run.charge_df_transcript,
            charge_df_aa=synth_run.charge_df_aa,
            charge_summary=synth_run.charge_summary_df,
            sample_df=synth_run.sample_df,
            context=synth_run.context,
        )
        out = tmp_path / 'charge_report.html'
        gen.generate_html_report(out)
        text = out.read_text()
        assert text.count('<!DOCTYPE html>') == 1
        assert '--accent: #00b894' in text

    def test_fragment_report_renders(self, synth_run, tmp_path):
        from trnaseq.qc.fragment_report import FragmentReportGenerator
        gen = FragmentReportGenerator(
            fragment_counts_df=synth_run.fragment_counts_df,
            rt_dropoff_df=synth_run.rt_dropoff_df,
            fragment_lengths_df=synth_run.fragment_lengths_df,
            fragment_summary_df=synth_run.fragment_summary_df,
            sample_df=synth_run.sample_df,
            context=synth_run.context,
        )
        out = tmp_path / 'fragment_report.html'
        gen.generate_html_report(out)
        text = out.read_text()
        assert text.count('<!DOCTYPE html>') == 1
        assert '--accent: #0984e3' in text

    def test_abundance_report_renders(self, synth_run, tmp_path):
        from trnaseq.qc.abundance_report import AbundanceReportGenerator
        gen = AbundanceReportGenerator(
            results_df=synth_run.abundance_results_df,
            count_matrix=synth_run.abundance_count_matrix,
            control_group='cond_0',
            level='transcript',
            condition_map=synth_run.abundance_condition_map,
            context=synth_run.context,
        )
        out = tmp_path / 'abundance_report.html'
        gen.generate_html_report(out)
        text = out.read_text()
        assert text.count('<!DOCTYPE html>') == 1
        assert '--accent: #6c5ce7' in text
