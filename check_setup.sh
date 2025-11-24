#!/bin/bash
# Quick test script to verify Nextflow pipeline setup

echo "============================================"
echo "Nextflow Pipeline Setup Verification"
echo "============================================"
echo ""

# Check if required files exist
echo "1. Checking required files..."
files=(
    "main.nf"
    "nextflow.config"
    "bin/run_annotation.py"
    "bin/lib/llm_client.py"
    "bin/lib/annotator.py"
    "bin/lib/data_processing.py"
    "bin/lib/logger.py"
    "data/test.h5ad"
)

all_good=true
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file (MISSING)"
        all_good=false
    fi
done

echo ""
echo "2. Checking Python script permissions..."
if [ -x "bin/run_annotation.py" ]; then
    echo "  ✓ bin/run_annotation.py is executable"
else
    echo "  ✗ bin/run_annotation.py is not executable"
    echo "    Run: chmod +x bin/run_annotation.py"
    all_good=false
fi

echo ""
echo "3. Checking Nextflow installation..."
if command -v nextflow &> /dev/null; then
    version=$(nextflow -version | head -1)
    echo "  ✓ Nextflow found: $version"
else
    echo "  ✗ Nextflow not found"
    echo "    Install from: https://www.nextflow.io/docs/latest/getstarted.html"
    all_good=false
fi

echo ""
echo "4. Checking Python environment..."
if command -v python &> /dev/null; then
    python_version=$(python --version 2>&1)
    echo "  ✓ Python found: $python_version"
    
    # Check for required packages
    echo "  Checking Python packages..."
    packages=("scanpy" "pandas" "numpy" "google.genai")
    for pkg in "${packages[@]}"; do
        if python -c "import $pkg" 2>/dev/null; then
            echo "    ✓ $pkg"
        else
            echo "    ✗ $pkg (MISSING)"
            all_good=false
        fi
    done
else
    echo "  ✗ Python not found"
    all_good=false
fi

echo ""
echo "============================================"
if [ "$all_good" = true ]; then
    echo "✅ Setup verification PASSED!"
    echo ""
    echo "You can now run the pipeline with:"
    echo "  nextflow run main.nf --input_h5ad data/test.h5ad --mock_llm -profile test"
else
    echo "⚠️  Setup verification FAILED"
    echo "Please fix the issues above before running the pipeline."
fi
echo "============================================"
