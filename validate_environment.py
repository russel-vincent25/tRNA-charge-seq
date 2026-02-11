#!/usr/bin/env python3
"""
Environment Validation Script for tRNA-charge-seq Pipeline

This script verifies that all required dependencies are correctly installed
and compatible. Run this after creating the conda environment.

Usage:
    conda activate tRNA-seq
    python validate_environment.py

Exit codes:
    0 - All checks passed
    1 - One or more checks failed
"""

import sys
import subprocess
from importlib.metadata import version, PackageNotFoundError

def check_python_package(package_name, min_version=None, max_version=None):
    """Check if a Python package is installed and meets version constraints."""
    try:
        installed_version = version(package_name)

        # Check minimum version
        if min_version and installed_version < min_version:
            print(f"❌ {package_name}: {installed_version} (need >={min_version})")
            return False

        # Check maximum version
        if max_version and installed_version >= max_version:
            print(f"❌ {package_name}: {installed_version} (need <{max_version})")
            return False

        print(f"✅ {package_name}: {installed_version}")
        return True

    except PackageNotFoundError:
        print(f"❌ {package_name}: NOT INSTALLED")
        return False

def check_system_tool(command, name=None):
    """Check if a system tool is available."""
    if name is None:
        name = command

    try:
        result = subprocess.run(
            [command, '--version'],
            capture_output=True,
            text=True,
            timeout=5
        )

        if result.returncode == 0:
            # Extract version from output (first line usually contains version)
            version_line = result.stdout.split('\n')[0] if result.stdout else result.stderr.split('\n')[0]
            print(f"✅ {name}: {version_line.strip()}")
            return True
        else:
            print(f"❌ {name}: Command failed")
            return False

    except FileNotFoundError:
        print(f"❌ {name}: NOT FOUND")
        return False
    except subprocess.TimeoutExpired:
        print(f"❌ {name}: Timeout")
        return False
    except Exception as e:
        print(f"❌ {name}: Error - {e}")
        return False

def check_imports():
    """Test critical imports that failed during beta testing."""
    print("\n=== Testing Critical Imports ===")

    try:
        import seaborn as sns
        print("✅ seaborn imports successfully (scipy compatibility OK)")
        success = True
    except Exception as e:
        print(f"❌ seaborn import failed: {e}")
        print("   → This usually means scipy version is incompatible")
        print("   → Try: pip install scipy==1.11.4 seaborn==0.13.0")
        success = False

    try:
        import pandas as pd
        df = pd.DataFrame({'a': [1, 2, 3]})
        print("✅ pandas imports and works")
    except Exception as e:
        print(f"❌ pandas import/usage failed: {e}")
        success = False

    return success

def main():
    print("=" * 70)
    print("tRNA-charge-seq Environment Validation")
    print("=" * 70)

    all_checks_passed = True

    # Python packages with version constraints from beta testing
    print("\n=== Python Packages ===")

    packages = [
        # Critical packages with version constraints (Issues #1 and #2)
        ('numpy', '1.24.0', '1.26.0'),
        ('pandas', '2.0.0', '2.2.0'),
        ('scipy', '1.11.0', '1.15.0'),
        ('seaborn', '0.13.0', None),

        # Other required packages
        ('matplotlib', '3.7.0', None),
        ('biopython', '1.80', None),
        ('jupyterlab', '4.0', None),
        ('logomaker', '0.8', None),
        ('mpire', '2.0', None),
        ('openpyxl', None, None),
        ('wand', None, None),
        ('natsort', None, None),
    ]

    for package_info in packages:
        package_name = package_info[0]
        min_ver = package_info[1] if len(package_info) > 1 else None
        max_ver = package_info[2] if len(package_info) > 2 else None

        if not check_python_package(package_name, min_ver, max_ver):
            all_checks_passed = False

    # System tools
    print("\n=== Bioinformatics Tools ===")

    tools = [
        ('AdapterRemoval', 'AdapterRemoval'),
        ('swipe', 'SWIPE'),
        ('makeblastdb', 'BLAST'),
        ('convert', 'ImageMagick'),
    ]

    for command, name in tools:
        if not check_system_tool(command, name):
            all_checks_passed = False

    # Critical import tests
    if not check_imports():
        all_checks_passed = False

    # Final summary
    print("\n" + "=" * 70)
    if all_checks_passed:
        print("✅ ALL CHECKS PASSED - Environment is ready!")
        print("=" * 70)
        print("\nYou can now run the tRNA-charge-seq pipeline.")
        print("Start with: jupyter lab projects/example/process_data_update.ipynb")
        return 0
    else:
        print("❌ SOME CHECKS FAILED - Please fix the issues above")
        print("=" * 70)
        print("\nRecommended fixes:")
        print("1. Recreate environment: conda env remove -n tRNA-seq")
        print("2. Install from environment.yml: conda env create -f environment.yml")
        print("3. Activate: conda activate tRNA-seq")
        print("4. Re-run this script: python validate_environment.py")
        return 1

if __name__ == "__main__":
    sys.exit(main())
