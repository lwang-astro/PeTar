version=`git rev-list --count HEAD`
version=`expr $version + 1`
echo $version >../VERSION
