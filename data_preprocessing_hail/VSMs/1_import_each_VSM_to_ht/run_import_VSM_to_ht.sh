#!/bin/bash

for f in import_*.py; do
  [ -e "$f" ] || continue
  python "$f" &
done