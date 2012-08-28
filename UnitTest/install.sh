#! /bin/sh
# Filename: install.sh
# You must set 'DOUBAN_LINALG_DIR' varaible before runing this file and notice there may be no '/' at the end.

if [ "$DOUBAN_LINALG_DIR" ]; then
  mv "$DOUBAN_LINALG_DIR/svd_tr.hpp" "$DOUBAN_LINALG_DIR/svd_tr.hpp.bak"
  mv "bidiag.hpp" "bidiag.hpp.bak"
  cp "svd_tr.hpp" "$DOUBAN_LINALG_DIR"
  cp "bidiag.hpp" "$DOUBAN_LINALG_DIR"
else
  echo "You must set 'DOUBAN_LINALG_DIR' varaible before running this script."
  exit 1
fi

echo "install svd test program finished!"
