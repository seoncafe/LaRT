#!/bin/bash
# Compile all .tex files in this directory and remove byproducts.
# Keeps: .tex, .cls, .pdf, .sh

DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$DIR"

for tex in *.tex; do
    echo "==> Compiling $tex ..."
    pdflatex -interaction=nonstopmode "$tex"
    pdflatex -interaction=nonstopmode "$tex"
    pdflatex -interaction=nonstopmode "$tex"
done

echo "==> Cleaning byproducts ..."
find . -maxdepth 1 -type f \
    ! -name "*.tex" \
    ! -name "*.cls" \
    ! -name "*.pdf" \
    ! -name "*.sh"  \
    -delete

echo "==> Done."
