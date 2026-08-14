#!/bin/bash
# Download all 18 tutorial SCA files via HTTPS

BASE_URL="https://stpubdata.s3.amazonaws.com/roman/nexus/soc_simulations/tutorial_data/roman-2026.2"
CACHE_DIR="cache"

mkdir -p "$CACHE_DIR"

SCAS=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18)

for sca in "${SCAS[@]}"; do
    FILENAME="r0003201001001001004_0001_wfi${sca}_f106_cal.asdf"
    LOCAL_PATH="$CACHE_DIR/$FILENAME"

    # Skip if already exists and is valid
    if [ -f "$LOCAL_PATH" ] && [ $(stat -f%z "$LOCAL_PATH") -gt 100000000 ]; then
        echo "✓ $FILENAME already cached ($(du -h "$LOCAL_PATH" | cut -f1))"
        continue
    fi

    URL="$BASE_URL/$FILENAME"
    echo "Downloading $FILENAME..."
    curl -o "$LOCAL_PATH" "$URL"

    if [ $? -eq 0 ]; then
        SIZE=$(du -h "$LOCAL_PATH" | cut -f1)
        echo "✓ Cached $FILENAME ($SIZE)"
    else
        echo "✗ Failed to download $FILENAME"
        rm -f "$LOCAL_PATH"
    fi
done

echo "Done. Cache directory contents:"
ls -lh "$CACHE_DIR" | tail -20
