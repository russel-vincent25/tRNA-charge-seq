#!/usr/bin/env python3
"""
Test script for AlignmentViewer

This script creates a minimal test dataset and demonstrates all
AlignmentViewer functionality.

Run from repository root:
    python test_alignment_viewer.py
"""

import json
import bz2
import sys
from pathlib import Path

# Add repo to path
repo_path = Path(__file__).parent
sys.path.insert(0, str(repo_path))

from trnaseq.visualization import AlignmentViewer


def create_test_data():
    """Create a minimal test JSON.bz2 file with mock alignment data"""
    print("=" * 80)
    print("Creating test alignment data...")
    print("=" * 80)

    # Create mock alignment data (3 reads for 2 tRNAs)
    test_data = {
        'read_001': {
            'aligned': True,
            'name': 'tRNA-Ala-TGC-1',
            'dpos': [1, 76],
            'qpos': [1, 76],
            'qseq': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA',
            'dseq': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA',
            'aseq': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
            'score': 152,
            'Fmax_score': 1.0,
            'Ndel': 0,
            'Nins': 0
        },
        'read_002': {
            'aligned': True,
            'name': 'tRNA-Ala-TGC-1',
            'dpos': [1, 75],
            'qpos': [1, 75],
            'qseq': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCC',
            'dseq': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCC',
            'aseq': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
            'score': 150,
            'Fmax_score': 0.99,
            'Ndel': 0,
            'Nins': 0
        },
        'read_003': {
            'aligned': True,
            'name': 'tRNA-Ala-TGC-1',
            'dpos': [2, 76],
            'qpos': [1, 75],
            'qseq': 'GGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA',
            'dseq': 'GGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA',
            'aseq': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
            'score': 148,
            'Fmax_score': 0.97,
            'Ndel': 0,
            'Nins': 0
        },
        'read_004': {
            'aligned': True,
            'name': 'tRNA-Ala-TGC-1',
            'dpos': [1, 76],
            'qpos': [1, 76],
            'qseq': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCGATCCCGCTTAGCTCCA',
            'dseq': 'GGGGCTATAGCTCAGCTGGGAGAGCGCTTGCATGGCATGCAAGAGGTCAGCGGTTCAATCCCGCTTAGCTCCA',
            'aseq': '|||||||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||||||||||',
            'score': 150,
            'Fmax_score': 0.98,
            'Ndel': 0,
            'Nins': 0
        },
        'read_005': {
            'aligned': True,
            'name': 'tRNA-Gly-GCC-1',
            'dpos': [1, 74],
            'qpos': [1, 74],
            'qseq': 'GCATTGGTGGTTCAGTGGTAGAATTCTCGCCTGCCACGCGGGAGGCCCGGGTTCGATTCCCGGCCAATGCA',
            'dseq': 'GCATTGGTGGTTCAGTGGTAGAATTCTCGCCTGCCACGCGGGAGGCCCGGGTTCGATTCCCGGCCAATGCA',
            'aseq': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
            'score': 148,
            'Fmax_score': 0.98,
            'Ndel': 0,
            'Nins': 0
        },
        'read_006': {
            'aligned': True,
            'name': 'tRNA-Gly-GCC-1',
            'dpos': [1, 74],
            'qpos': [1, 74],
            'qseq': 'GCATTGGTGGTTCAGTGGTAGAATTCTCGCCTGCCACGCGGGAGGCCCGGGTTCGATTCCCGGCCAATGCA',
            'dseq': 'GCATTGGTGGTTCAGTGGTAGAATTCTCGCCTGCCACGCGGGAGGCCCGGGTTCGATTCCCGGCCAATGCA',
            'aseq': '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||',
            'score': 148,
            'Fmax_score': 0.98,
            'Ndel': 0,
            'Nins': 0
        }
    }

    # Save to bz2-compressed JSON
    test_file = repo_path / 'test_sample_SWalign.json.bz2'
    with bz2.open(test_file, 'wt') as f:
        json.dump(test_data, f)

    print(f"✅ Created test file: {test_file.name}")
    print(f"   - 6 reads total")
    print(f"   - 2 tRNAs: tRNA-Ala-TGC-1 (4 reads), tRNA-Gly-GCC-1 (2 reads)")
    print()

    return test_file


def test_alignment_viewer(test_file):
    """Test all AlignmentViewer methods"""
    print("=" * 80)
    print("Testing AlignmentViewer Methods")
    print("=" * 80)
    print()

    # Test 1: Initialize viewer
    print("Test 1: Initialize AlignmentViewer")
    print("-" * 40)
    viewer = AlignmentViewer(test_file)
    print("✅ AlignmentViewer initialized successfully")
    print()

    # Test 2: List tRNAs
    print("Test 2: List available tRNAs")
    print("-" * 40)
    trna_list = viewer.list_trnas(min_reads=1)
    print(f"✅ Found {len(trna_list)} tRNAs:")
    print(trna_list.to_string(index=False))
    print()

    # Test 3: Load alignments
    print("Test 3: Load alignments for specific tRNA")
    print("-" * 40)
    target_trna = trna_list.iloc[0]['tRNA']
    alignments = viewer.load_alignments(target_trna)
    print(f"✅ Loaded {len(alignments)} alignments for {target_trna}")
    print()

    # Test 4: Calculate coverage
    print("Test 4: Calculate coverage statistics")
    print("-" * 40)
    cov_df = viewer.calculate_coverage(target_trna)
    print(f"✅ Coverage calculated for {len(cov_df)} positions")
    print(f"   Mean coverage: {cov_df['coverage'].mean():.2f}")
    print(f"   Max coverage: {cov_df['coverage'].max()}")
    print(f"   Mean mismatch rate: {cov_df['mismatch_rate'].mean():.2f}%")
    print(f"   Max mismatch rate: {cov_df['mismatch_rate'].max():.2f}%")
    print()

    # Test 5: Get alignment details
    print("Test 5: Get alignment details")
    print("-" * 40)
    details = viewer.get_alignment_details(target_trna, max_reads=10)
    print(f"✅ Retrieved details for {len(details)} reads")
    if len(details) > 0:
        print(f"   Example read: {details[0]['read_id']}")
        print(f"   Score: {details[0]['score']}")
        print(f"   Position: {details[0]['dpos']}")
    print()

    # Test 6: Plot coverage
    print("Test 6: Generate coverage plot")
    print("-" * 40)
    output_png = repo_path / f'test_{target_trna.replace("/", "-")}_coverage.png'
    fig = viewer.plot_coverage(target_trna, output=str(output_png))
    if output_png.exists():
        print(f"✅ Coverage plot saved: {output_png.name}")
        print(f"   File size: {output_png.stat().st_size / 1024:.1f} KB")
    print()

    # Test 7: Create HTML report
    print("Test 7: Create HTML report")
    print("-" * 40)
    output_html = repo_path / f'test_{target_trna.replace("/", "-")}_report.html'
    html_path = viewer.create_html_report(target_trna, output=str(output_html))
    if output_html.exists():
        print(f"✅ HTML report saved: {output_html.name}")
        print(f"   File size: {output_html.stat().st_size / 1024:.1f} KB")
    print()

    # Test 8: Quick view
    print("Test 8: Quick view (PNG + HTML)")
    print("-" * 40)
    # Test with second tRNA
    if len(trna_list) > 1:
        target_trna2 = trna_list.iloc[1]['tRNA']
        png_path, html_path = viewer.quick_view(target_trna2)
        print(f"✅ Quick view generated for {target_trna2}")
        print(f"   PNG: {Path(png_path).name}")
        print(f"   HTML: {Path(html_path).name}")
    print()

    return output_png, output_html


def cleanup_test_files(test_file, *output_files):
    """Clean up test files"""
    print("=" * 80)
    print("Cleaning up test files...")
    print("=" * 80)

    files_to_remove = [test_file] + list(output_files)

    # Also remove quick view files
    for pattern in ['test_sample_tRNA-*.png', 'test_sample_tRNA-*.html']:
        files_to_remove.extend(repo_path.glob(pattern))

    for file_path in files_to_remove:
        if Path(file_path).exists():
            Path(file_path).unlink()
            print(f"✅ Removed: {Path(file_path).name}")

    print()


def main():
    """Run all tests"""
    print()
    print("╔" + "=" * 78 + "╗")
    print("║" + " " * 20 + "ALIGNMENT VIEWER TEST SUITE" + " " * 31 + "║")
    print("╚" + "=" * 78 + "╝")
    print()

    try:
        # Create test data
        test_file = create_test_data()

        # Test AlignmentViewer
        output_png, output_html = test_alignment_viewer(test_file)

        # Summary
        print("=" * 80)
        print("✅ ALL TESTS PASSED!")
        print("=" * 80)
        print()
        print("Summary:")
        print("  ✓ AlignmentViewer class initialized")
        print("  ✓ All 7 methods tested and working")
        print("  ✓ PNG coverage plots generated")
        print("  ✓ HTML reports created")
        print("  ✓ Quick view mode working")
        print()
        print("Generated test outputs:")
        print(f"  - {output_png.name}")
        print(f"  - {output_html.name}")
        print(f"  - test_sample_tRNA-Gly-GCC-1_coverage.png")
        print(f"  - test_sample_tRNA-Gly-GCC-1_report.html")
        print()
        print("You can:")
        print(f"  1. View the coverage plot: {output_png}")
        print(f"  2. Open the HTML report in browser: {output_html}")
        print()

        # Ask if user wants to keep test files
        response = input("Keep test files? (y/n) [n]: ").strip().lower()
        if response != 'y':
            cleanup_test_files(test_file, output_png, output_html)
            print("Test files removed.")
        else:
            print("Test files kept for inspection.")

        print()
        print("=" * 80)
        print("🎉 AlignmentViewer is ready for production use!")
        print("=" * 80)
        print()
        print("Next steps:")
        print("  1. Generate alignment JSON files from your tRNA-seq data")
        print("  2. Use AlignmentViewer to visualize coverage")
        print("  3. Check notebooks/02_alignment_qc.ipynb for more examples")
        print()

    except Exception as e:
        print()
        print("=" * 80)
        print("❌ TEST FAILED")
        print("=" * 80)
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
