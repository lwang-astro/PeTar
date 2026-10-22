version=`git rev-list --count HEAD`
version=`expr $version + 1`
# Preserve the 'e' experiment suffix across bumps (same fix as SDAR's
# get_version.sh, 2026-10-22: the manual append step was error-prone).
suffix=''
[ -f ../VERSION ] && grep -q 'e$' ../VERSION && suffix='e'
echo "${version}${suffix}" >../VERSION
