#!/bin/bash
# Starts Geneious to GFF3 version 1.1. Double-click this file in Finder.
# Nothing needs to be installed except Python: Flask is included in the "vendor" folder.
cd "$(dirname "$0")" || exit 1
# On a Mac the jobs live in this folder, so version 1.0 in its own folder is untouched.
: "${GFF3TOOL_PORT:=8775}"
export GFF3TOOL_PORT
if command -v python3 >/dev/null 2>&1 && python3 -c "import sys; sys.exit(sys.version_info < (3, 9))" 2>/dev/null; then
  exec python3 app.py
fi
echo "Python 3.9 or newer was not found. Install it from https://www.python.org/downloads/ and try again."
read -r -p "Press Enter to close."
exit 1
