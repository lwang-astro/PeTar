./configure 
make clean
make install
./configure --with-interrupt=bse
make clean
make install
./configure --with-interrupt=bse --with-external=galpy
make clean
make install
./configure --with-external=galpy --enable-mpfrc
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
./configure --with-interrupt=bse --with-external=galpy --with-external-hard=gasdrag
make clean
make install
./configure --with-interrupt=bse --with-external-hard=gasdrag
make clean
make install
./configure --with-interrupt=dsm --with-external=galpy --with-external-hard=gasdrag
make clean
make install
./configure --with-interrupt=dsm --with-external-hard=gasdrag
make clean
make install
./configure --with-external-hard=gasdrag
make clean
make install
