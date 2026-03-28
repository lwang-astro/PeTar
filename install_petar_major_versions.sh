./configure 
make clean
make install
./configure --with-interrupt=bse
make clean
make install
./configure --with-interrupt=bse --with-external=galpy
make clean
make install
./configure --with-external=galpy
make clean
make install
./configure --with-interrupt=bse --with-external=agama
make clean
make install
./configure --with-external=agama
make clean
make install
