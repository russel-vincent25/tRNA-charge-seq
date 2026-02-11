"""
Lightweight Alignment Viewer for tRNA-seq JSON Data

This module provides a lightweight alignment viewer that reads JSON.bz2 files
directly and generates IGV-like visualizations without needing external tools.

Classes:
    AlignmentViewer: Main class for alignment visualization and reporting

Usage:
    from trnaseq.visualization import AlignmentViewer

    viewer = AlignmentViewer('sample001_SWalign.json.bz2')
    viewer.plot_coverage('tRNA-Ala-TGC-1', output='ala_coverage.png')
    viewer.create_html_report('tRNA-Ala-TGC-1', output='ala_report.html')
"""

import bz2
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import base64
from io import BytesIO


class AlignmentViewer:
    """
    Lightweight alignment viewer for tRNA-seq JSON data

    Creates IGV-like visualizations directly from JSON.bz2 files
    without needing external tools.

    Usage:
        viewer = AlignmentViewer('sample001_SWalign.json.bz2')
        viewer.plot_coverage('tRNA-Ala-TGC-1', output='ala_coverage.png')
        viewer.create_html_report('tRNA-Ala-TGC-1', output='ala_report.html')
    """

    def __init__(self, json_path, reference_seqs=None):
        """
        Initialize the AlignmentViewer

        Parameters:
            json_path: Path to {sample}_SWalign.json.bz2 file
            reference_seqs: Dict of {trna_id: sequence} (optional)
        """
        self.json_path = Path(json_path)
        self.reference_seqs = reference_seqs or {}
        self.alignments = None

        if not self.json_path.exists():
            raise FileNotFoundError(f"JSON file not found: {self.json_path}")

    def load_alignments(self, trna_id=None):
        """
        Load alignments from JSON file

        Parameters:
            trna_id: Load only alignments for this tRNA (memory efficient)
                     If None, load all alignments

        Returns:
            Dictionary of alignments {read_id: alignment_data}
        """
        print(f"Loading alignments from {self.json_path.name}...")

        with bz2.open(self.json_path, 'rt') as f:
            all_alignments = json.load(f)

        if trna_id:
            # Filter for specific tRNA
            filtered = {
                read_id: aln
                for read_id, aln in all_alignments.items()
                if aln.get('aligned') and aln.get('name') == trna_id
            }
            print(f"  Loaded {len(filtered)} reads aligned to {trna_id}")
            return filtered
        else:
            print(f"  Loaded {len(all_alignments)} total reads")
            return all_alignments

    def calculate_coverage(self, trna_id, ref_length=None):
        """
        Calculate per-position coverage and mismatch rates

        Parameters:
            trna_id: tRNA annotation to analyze
            ref_length: Expected reference length (auto-detected if None)

        Returns:
            DataFrame with columns: position, coverage, mismatches, mismatch_rate
        """
        alignments = self.load_alignments(trna_id)

        if not alignments:
            print(f"No alignments found for {trna_id}")
            return None

        # Determine reference length
        if ref_length is None:
            # Use max alignment end position
            ref_length = max(
                aln['dpos'][1]
                for aln in alignments.values()
                if aln.get('aligned')
            )

        # Initialize coverage arrays
        coverage = np.zeros(ref_length + 1, dtype=int)
        mismatches = np.zeros(ref_length + 1, dtype=int)

        # Count coverage and mismatches
        for aln in alignments.values():
            if not aln.get('aligned'):
                continue

            # Alignment positions on reference (1-indexed)
            start, end = aln['dpos']

            # Increment coverage
            coverage[start:end+1] += 1

            # Count mismatches from alignment strings
            qseq = aln['qseq']  # Query sequence
            dseq = aln['dseq']  # Database (reference) sequence

            # Track position on reference
            ref_pos = start
            for i, (q, d) in enumerate(zip(qseq, dseq)):
                # Skip gaps in reference
                if d == '-':
                    continue

                # Count mismatch (not gap, and not match)
                if q != d and q != '-':
                    if ref_pos <= ref_length:
                        mismatches[ref_pos] += 1

                ref_pos += 1

        # Create DataFrame
        df = pd.DataFrame({
            'position': range(1, ref_length + 1),
            'coverage': coverage[1:],  # Skip position 0
            'mismatches': mismatches[1:],
            'mismatch_rate': np.divide(
                mismatches[1:],
                coverage[1:],
                out=np.zeros(ref_length),
                where=coverage[1:] > 0
            ) * 100
        })

        return df

    def plot_coverage(self, trna_id, output=None, figsize=(14, 6)):
        """
        Create IGV-like coverage plot

        Parameters:
            trna_id: tRNA annotation to plot
            output: Output file path (PNG, PDF, or None for display)
            figsize: Figure size (width, height)

        Returns:
            matplotlib Figure object
        """
        cov_df = self.calculate_coverage(trna_id)

        if cov_df is None:
            return None

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize,
                                        height_ratios=[3, 1])

        # Top panel: Coverage with mismatches
        ax1.bar(cov_df['position'], cov_df['coverage'],
                color='gray', alpha=0.7, label='Coverage')
        ax1.bar(cov_df['position'], cov_df['mismatches'],
                color='red', alpha=0.8, label='Mismatches')
        ax1.set_ylabel('Read Count', fontsize=12)
        ax1.set_title(f'Alignment Coverage: {trna_id}', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(axis='y', alpha=0.3)

        # Bottom panel: Mismatch rate heatmap
        sns.heatmap(
            cov_df['mismatch_rate'].values.reshape(1, -1),
            cmap='Reds', cbar_kws={'label': 'Mismatch %'},
            xticklabels=False, yticklabels=False, ax=ax2
        )
        ax2.set_xlabel('Position on tRNA', fontsize=12)
        ax2.set_ylabel('Mismatch\nRate', fontsize=10)

        plt.tight_layout()

        if output:
            plt.savefig(output, dpi=300, bbox_inches='tight')
            print(f"  Saved: {output}")
        else:
            plt.show()

        return fig

    def get_alignment_details(self, trna_id, max_reads=50):
        """
        Get detailed alignment information for visualization

        Parameters:
            trna_id: tRNA annotation
            max_reads: Maximum number of reads to return (for display)

        Returns:
            List of dicts with alignment details
        """
        alignments = self.load_alignments(trna_id)

        details = []
        for i, (read_id, aln) in enumerate(alignments.items()):
            if i >= max_reads:
                break

            if not aln.get('aligned'):
                continue

            details.append({
                'read_id': read_id,
                'score': aln['score'],
                'dpos': aln['dpos'],
                'qpos': aln.get('qpos', [0, 0]),
                'qseq': aln['qseq'],
                'aseq': aln['aseq'],  # Alignment string (| for match)
                'dseq': aln['dseq'],
                'Ndel': aln.get('Ndel', 0),
                'Nins': aln.get('Nins', 0),
                'Fmax_score': aln.get('Fmax_score', 0)
            })

        return details

    def create_html_report(self, trna_id, output='alignment_report.html'):
        """
        Create interactive HTML report with coverage plot and read details

        Parameters:
            trna_id: tRNA annotation
            output: Output HTML file path

        Returns:
            Path to output HTML file
        """
        # Generate coverage plot
        fig = self.plot_coverage(trna_id, output=None)

        if fig is None:
            print(f"Cannot create report: No alignments found for {trna_id}")
            return None

        # Convert plot to base64 for embedding in HTML
        buf = BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        plt.close(fig)

        # Get alignment details
        details = self.get_alignment_details(trna_id, max_reads=100)

        # Generate HTML
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Alignment Report: {trna_id}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #fafafa; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .coverage-plot {{ max-width: 100%; margin: 20px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .alignment-table {{
            border-collapse: collapse;
            width: 100%;
            margin-top: 20px;
            font-family: monospace;
            font-size: 12px;
            background: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .alignment-table th {{
            background-color: #3498db;
            color: white;
            padding: 10px;
            text-align: left;
            position: sticky;
            top: 0;
        }}
        .alignment-table td {{
            padding: 8px;
            border-bottom: 1px solid #ddd;
        }}
        .alignment-table tr:hover {{ background-color: #f5f5f5; }}
        .alignment-seq {{
            font-family: 'Courier New', monospace;
            white-space: pre;
            background: #f8f9fa;
            padding: 5px;
            border-radius: 3px;
            font-size: 11px;
        }}
        .stats {{
            background: white;
            padding: 15px;
            border-radius: 5px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 10px;
        }}
        .stat-item {{
            padding: 10px;
            background: #ecf0f1;
            border-radius: 3px;
        }}
        .stat-label {{
            font-weight: bold;
            color: #7f8c8d;
            font-size: 12px;
        }}
        .stat-value {{
            font-size: 18px;
            color: #2c3e50;
            margin-top: 5px;
        }}
    </style>
</head>
<body>
    <h1>Alignment Report: {trna_id}</h1>

    <div class="stats">
        <div class="stats-grid">
            <div class="stat-item">
                <div class="stat-label">Sample</div>
                <div class="stat-value">{self.json_path.stem.replace('_SWalign', '')}</div>
            </div>
            <div class="stat-item">
                <div class="stat-label">tRNA</div>
                <div class="stat-value">{trna_id}</div>
            </div>
            <div class="stat-item">
                <div class="stat-label">Aligned Reads</div>
                <div class="stat-value">{len(details)}</div>
            </div>
        </div>
    </div>

    <h2>Coverage Plot</h2>
    <img src="data:image/png;base64,{img_base64}" class="coverage-plot">

    <h2>Read Alignments (Top {len(details)} reads)</h2>
    <table class="alignment-table">
        <tr>
            <th>Read ID</th>
            <th>Score</th>
            <th>Ref Pos</th>
            <th>Alignment</th>
            <th>Fmax</th>
        </tr>
"""

        for aln in details:
            # Format alignment strings for display
            alignment_display = f"""Query: {aln['qseq']}
       {aln['aseq']}
Ref:   {aln['dseq']}"""

            # Truncate read ID for display
            read_id_display = aln['read_id'][:30] + '...' if len(aln['read_id']) > 30 else aln['read_id']

            html += f"""        <tr>
            <td>{read_id_display}</td>
            <td>{aln['score']}</td>
            <td>{aln['dpos'][0]}-{aln['dpos'][1]}</td>
            <td class="alignment-seq">{alignment_display}</td>
            <td>{aln['Fmax_score']:.2f}</td>
        </tr>
"""

        html += """    </table>
</body>
</html>
"""

        # Write HTML file
        with open(output, 'w') as f:
            f.write(html)

        print(f"\nHTML report created: {output}")
        print(f"   Open in browser to view alignments")

        return output

    def quick_view(self, trna_id):
        """
        Quick view: Generate both PNG and HTML in one command

        Parameters:
            trna_id: tRNA annotation to view

        Returns:
            Tuple of (png_path, html_path)
        """
        sample_name = self.json_path.stem.replace('_SWalign', '')

        # Generate PNG
        png_output = f"{sample_name}_{trna_id}_coverage.png"
        self.plot_coverage(trna_id, output=png_output)

        # Generate HTML report
        html_output = f"{sample_name}_{trna_id}_report.html"
        self.create_html_report(trna_id, output=html_output)

        print(f"\nQuick view complete!")
        print(f"   PNG:  {png_output}")
        print(f"   HTML: {html_output}")

        return png_output, html_output

    def list_trnas(self, min_reads=10):
        """
        List all tRNAs present in the alignment file with read counts

        Parameters:
            min_reads: Minimum number of reads to include in listing

        Returns:
            DataFrame with tRNA names and read counts
        """
        print(f"Scanning {self.json_path.name} for tRNA annotations...")

        with bz2.open(self.json_path, 'rt') as f:
            all_alignments = json.load(f)

        # Count reads per tRNA
        trna_counts = {}
        for read_id, aln in all_alignments.items():
            if aln.get('aligned'):
                trna_name = aln.get('name')
                if trna_name:
                    trna_counts[trna_name] = trna_counts.get(trna_name, 0) + 1

        # Create DataFrame
        df = pd.DataFrame([
            {'tRNA': name, 'reads': count}
            for name, count in trna_counts.items()
            if count >= min_reads
        ])

        # Sort by read count
        df = df.sort_values('reads', ascending=False).reset_index(drop=True)

        print(f"  Found {len(df)} tRNAs with >= {min_reads} reads")
        print(f"  Total tRNAs: {len(trna_counts)}")

        return df
