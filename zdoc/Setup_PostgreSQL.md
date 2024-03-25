# -------------------------------------------------------------------
# Create PostgreSQL with RDKit - Linux [imb-coadd-work]
# -------------------------------------------------------------------
# Based on https://iwatobipen.wordpress.com/2023/12/23/build-postgresql-rdkit-cartridge-and-install-new-version-of-postgresql-rdkit-postgresql-cheminformatics/
#

# Create Conda enviroment
conda create -n pg16 python
conda activate pg16

conda install -c conda-forge -c anaconda -y cmake cairo pillow eigen pkg-config
conda install -c conda-forge -c anaconda  -y boost-cpp boost py-boost
conda install -c conda-forge matplotlib gxx_linux-64 postgresql

conda list
# packages in environment at /opt/conda/envs/pg16:
#
# Name                    Version                   Build  Channel
binutils_linux-64         2.40                 hdade7a5_3    conda-forge
boost                     1.82.0          py312h06a4308_2    anaconda
cairo                     1.18.0               h3faef2a_0    conda-forge
cmake                     3.28.4               hcfe8598_0    conda-forge
eigen                     3.4.0                h00ab1b0_0    conda-forge
glib                      2.80.0               hf2295e7_1    conda-forge
gxx_linux-64              13.2.0               he8deefe_3    conda-forge
kernel-headers_linux-64   2.6.32              he073ed8_17    conda-forge
libblas                   3.9.0           21_linux64_openblas    conda-forge
matplotlib-base           3.8.3           py312he5832f3_0    conda-forge
openssl                   3.2.1                hd590300_1    conda-forge
postgresql                16.2                 h82ecc9d_1    conda-forge
python                    3.12.2          hab00c5b_0_cpython    conda-forge

# Compile rdkit-postgresql
cmake -DPy_ENABLE_SHARED=1 \
  -DRDK_INSTALL_INTREE=ON \
  -DRDK_INSTALL_STATIC_LIBS=OFF \
  -DRDK_BUILD_CPP_TESTS=ON \
  -DRDK_BUILD_PGSQL=ON \
  -DRDK_PGSQL_STATIC=ON \
  -DPYTHON_NUMPY_INCLUDE_PATH="$(python -c 'import numpy ; print(numpy.get_include())')" \
  -DBOOST_ROOT="$CONDA_PREFIX" \
  -DPostgreSQL_INCLUDE_DIR=$CONDA_PREFIX/include \
  -DPostgreSQL_TYPE_INCLUDE_DIR=$CONDA_PREFIX/include/server \
  -DPostgreSQL_LIBRARY=$CONDA_PREFIX/lib/libpq.so \
  ..
  
make
make install

# Test rdkit-postgresql
conda install pytest pandas
RDBASE=$PWD/.. PYTHONPATH=$RDBASE LD_LIBRARY_PATH=$RDBASE/lib:$LD_LIBRARY_PATH ctest

See also: https://github.com/rdkit/rdkit/discussions/6148

# -------------------------------------------------------------------
# Create First Database
# -------------------------------------------------------------------
initdb -D /opt/postgres/data/db_pg16 -U postgres
pg_ctl -D /opt/postgres/data/db_pg16 -l logfile start
psql -U postgres

# pg_ctl -D /opt/postgres/data/db_pg16 stop

# -------------------------------------------------------------------
# Create orgdb Tablespace/Database/User/Schema
# either by psql or pgAdmin
# -------------------------------------------------------------------

# Create orgdb USER
CREATE ROLE orgdb WITH
    LOGIN
    NOSUPERUSER
    NOCREATEDB
    NOCREATEROLE
    INHERIT
    NOREPLICATION
    CONNECTION LIMIT -1
    PASSWORD 'orgdb';

COMMENT ON ROLE orgdb IS 'orgdb';

# Create new TABLESPACE
mkdir /opt/postgres/data/<newDB>
chmod 700 /opt/postgres/data/<newDB>
CREATE TABLESPACE <newTbl>
  OWNER <newUser>
  LOCATION '/opt/postgres/data/<newDB>';

ALTER TABLESPACE <newTbl>
  OWNER TO <newUser>;
Database

# Create Database
CREATE DATABASE <newDB>
    WITH 
    OWNER = <newUser>
    ENCODING = 'UTF8'
    TABLESPACE = <newTbl>
    CONNECTION LIMIT = -1;
Add Extension

# Add RDKit Extension
CREATE EXTENSION rdkit;
Schema

# Create schemas
CREATE SCHEMA <newSchema>
    AUTHORIZATION <newDB>;
