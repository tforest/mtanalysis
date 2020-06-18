# ./write_taxonomy.sh PATH_PYTHON_SCRIPT DIR_WITH_GB_FILES
for f in $2/*; do
  python3 $1 "$f"
done
