"""
Example usage of AlignmentViewer

This script demonstrates basic usage of the AlignmentViewer class
for quick alignment quality control.

Run this script from the repository root:
    python -m trnaseq.visualization.example_usage
"""

import sys
from pathlib import Path

# Add repo to path if needed
repo_path = Path(__file__).parent.parent.parent
sys.path.insert(0, str(repo_path))

from trnaseq.visualization import AlignmentViewer


def example_basic_usage():
    """Example 1: Basic usage - load and visualize a single tRNA"""
    print("=" * 80)
    print("Example 1: Basic Usage")
    print("=" * 80)

    # NOTE: Replace this with your actual alignment file path
    json_file = repo_path / 'projects' / 'example' / 'data' / 'SWalign' / '100p_SWalign.json.bz2'

    if not json_file.exists():
        print(f"\nAlignment file not found: {json_file}")
        print("Please adjust the path to point to your alignment JSON file.")
        print("\nExpected structure:")
        print("  projects/example/data/SWalign/sample_SWalign.json.bz2")
        return False

    # Create viewer
    viewer = AlignmentViewer(json_file)

    # List available tRNAs
    print(f"\nListing tRNAs in {json_file.name}...")
    trna_list = viewer.list_trnas(min_reads=10)

    if len(trna_list) == 0:
        print("No tRNAs found with sufficient reads")
        return False

    print(f"\nTop 5 tRNAs by read count:")
    print(trna_list.head())

    # Visualize first tRNA
    selected_trna = trna_list.iloc[0]['tRNA']
    print(f"\nGenerating coverage plot for: {selected_trna}")

    # Plot coverage (saves to file)
    viewer.plot_coverage(selected_trna, output=f'example_{selected_trna.replace("/", "-")}_coverage.png')

    print(f"\nBasic usage complete!")
    return True


def example_quick_view():
    """Example 2: Quick view - generate both PNG and HTML"""
    print("\n" + "=" * 80)
    print("Example 2: Quick View (PNG + HTML)")
    print("=" * 80)

    json_file = repo_path / 'projects' / 'example' / 'data' / 'SWalign' / '100p_SWalign.json.bz2'

    if not json_file.exists():
        print(f"\nAlignment file not found: {json_file}")
        return False

    viewer = AlignmentViewer(json_file)
    trna_list = viewer.list_trnas(min_reads=10)

    if len(trna_list) == 0:
        print("No tRNAs found")
        return False

    selected_trna = trna_list.iloc[0]['tRNA']
    print(f"\nGenerating quick view for: {selected_trna}")

    # Quick view generates both PNG and HTML
    png_path, html_path = viewer.quick_view(selected_trna)

    print(f"\nQuick view complete!")
    print(f"  PNG:  {png_path}")
    print(f"  HTML: {html_path}")
    print(f"\nOpen {html_path} in your browser to view the interactive report!")
    return True


def example_coverage_analysis():
    """Example 3: Analyze coverage statistics"""
    print("\n" + "=" * 80)
    print("Example 3: Coverage Analysis")
    print("=" * 80)

    json_file = repo_path / 'projects' / 'example' / 'data' / 'SWalign' / '100p_SWalign.json.bz2'

    if not json_file.exists():
        print(f"\nAlignment file not found: {json_file}")
        return False

    viewer = AlignmentViewer(json_file)
    trna_list = viewer.list_trnas(min_reads=10)

    if len(trna_list) == 0:
        print("No tRNAs found")
        return False

    # Analyze top 3 tRNAs
    print(f"\nAnalyzing coverage for top 3 tRNAs...")

    for idx, row in trna_list.head(3).iterrows():
        trna_name = row['tRNA']
        read_count = row['reads']

        cov_df = viewer.calculate_coverage(trna_name)

        print(f"\n{idx+1}. {trna_name} (n={read_count} reads)")
        print(f"   Length: {len(cov_df)} positions")
        print(f"   Mean coverage: {cov_df['coverage'].mean():.1f}")
        print(f"   Max coverage: {cov_df['coverage'].max()}")
        print(f"   Mean mismatch rate: {cov_df['mismatch_rate'].mean():.2f}%")

        # Find high mismatch positions
        high_mm = cov_df[cov_df['mismatch_rate'] > 5.0]
        if len(high_mm) > 0:
            print(f"   Positions with >5% mismatch: {len(high_mm)}")
        else:
            print(f"   No positions with >5% mismatch (good quality!)")

    return True


def example_batch_processing():
    """Example 4: Batch process multiple samples"""
    print("\n" + "=" * 80)
    print("Example 4: Batch Processing")
    print("=" * 80)

    align_dir = repo_path / 'projects' / 'example' / 'data' / 'SWalign'

    if not align_dir.exists():
        print(f"\nAlignment directory not found: {align_dir}")
        return False

    # Find all alignment files
    json_files = list(align_dir.glob('*_SWalign.json.bz2'))

    if len(json_files) == 0:
        print(f"No alignment files found in {align_dir}")
        return False

    print(f"\nFound {len(json_files)} alignment files")

    # Process first 3 files as example
    target_trna = None

    for json_file in json_files[:3]:
        print(f"\nProcessing: {json_file.name}")

        viewer = AlignmentViewer(json_file)
        trna_list = viewer.list_trnas(min_reads=10)

        if len(trna_list) > 0:
            # Use first tRNA from first file as target
            if target_trna is None:
                target_trna = trna_list.iloc[0]['tRNA']

            # Check if target tRNA exists in this sample
            if target_trna in trna_list['tRNA'].values:
                sample_name = json_file.stem.replace('_SWalign', '')
                output_file = f"batch_{sample_name}_{target_trna.replace('/', '-')}.png"

                viewer.plot_coverage(target_trna, output=output_file)
                print(f"  Generated: {output_file}")
            else:
                print(f"  {target_trna} not found in this sample")

    print(f"\nBatch processing complete!")
    return True


def main():
    """Run all examples"""
    print("AlignmentViewer - Example Usage")
    print("=" * 80)

    examples = [
        example_basic_usage,
        example_quick_view,
        example_coverage_analysis,
        example_batch_processing
    ]

    success_count = 0
    for example_func in examples:
        try:
            if example_func():
                success_count += 1
        except Exception as e:
            print(f"\nError in {example_func.__name__}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 80)
    print(f"Examples complete: {success_count}/{len(examples)} successful")
    print("=" * 80)

    if success_count == 0:
        print("\nNote: Examples require alignment files to be present.")
        print("If you don't have alignment files yet, run the main processing")
        print("pipeline first to generate alignment JSON files.")
        print("\nSee: projects/example/process_data.ipynb")


if __name__ == '__main__':
    main()
