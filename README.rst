Codes to generate star catalogs for Stellarium
======================================================

This repository contains codes to generate experimental star catalogs for the popular open source planetarium software `Stellarium`_.
This experimental star catalogs are not meant to be used for astrophysical research, but rather aiming to provide a realistic 
sky for the Stellarium planetarium software.

The star catalogs are generated from the Hipparcos and Gaia DR3 catalogs.

.. _Stellarium: https://stellarium.org/

Requirements
------------

The code here are primarily written in Python using `Python>=3.12`. Run the following command to install the required Python packages:

..  code-block:: bash

   pip install -r requirements.txt


You also need to install `embree`_ for high performance ray tracing to quickly check which zones stars fall within on the sky.

.. _embree: https://www.embree.org/

Getting Data
------------

Most of the data are obtained from SIMBAD for Hipparcos stars and Gaia Archive for Gaia DR3 stars. To get Hipparcos stars from SIMBAD, run the following python script

.. code-block:: bash

   python simbad_query_hipsaohdhr.py.py

This will query SIMBAD for all possible Hipparcos IDs and save as a astropy Table.

To generate synthetic photometry (mainly `B-V` color for visual purposes) from Gaia DR3 continuous spectra, you need to download a large amount of data (~15 TB free disk space). The procedure is as follows:

#. | Download the data (DR3: 3362 GB) from this link - http://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/
   | To do this programmatically, you can run ``wget -P ./ --no-clobber --no-verbose --no-parent --recursive --level=1 --no-directories http://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/``
#. | Unzip the data (DR3: 8204 GB)
   | To do this programmatically with 8 CPU cores, you can run ``ls *.gz | xargs -n 1 -P 8 -I {} 'gunzip -c {} > ./$(basename {} .gz)'``
#. | Run GaiaXPy to synethize JKC photometry from the continuous spectra
   | To do this programmatically with 8 CPU cores, you can run ``ls ./XpContinuousMeanSpectrum_*.csv | xargs -n 1 -P 8 -I {} sh -c 'python gaiaxpy_phot.py {} ./jkc_photometry/tables'``
#. | After the above is completed and to merge into a single table with only `source_id`, $B$ and $V$ photometry, run ``python gaiaxpy_phot_merge.py``

To do ADQL query for Gaia DR3 stars, you can simply go to https://gea.esac.esa.int/archive/ and click `Search` and then `Advanced (ADQL)` and paste the query which should be inside Jupyter Notebooks in this repository.

Scripts and Notebooks
----------------------

Here are the scripts and notebooks in this repository, you should **follow the order** to run the scripts and notebooks.

#. | `1_simbad_query_hipsaohdhr.py`_ - This script query SIMBAD for all possible Hipparcos/SAO/HD/HR IDs and save as a astropy table.
#. | `Parse_HIP_Catalog.ipynb`_ - This notebook parse the HIP catalog and make sure the data is correct and clean along with Gaia source IDs and binary component IDs.
#. | `Gaia_Photometry.ipynb`_ - (*Optional*) This notebook provides sanity check for $V$-band magnitude and $B-V$ color computed from the Gaia DR3 photometry.
#. | `Star_Catalog_lv0_6.ipynb`_ - This notebook generate a bright star catalog (0-6) for Stellarium.
#. | `Star_catalog_lv7.ipynb`_ - This notebook generate a star catalog level 7 for Stellarium.
#. | `Star_catalog_lv8.ipynb`_ - This notebook generate a star catalog level 8 for Stellarium.

.. _1_simbad_query_hipsaohdhr.py: 1_simbad_query_hipsaohdhr.py
.. _Parse_HIP_Catalog.ipynb: Parse_HIP_Catalog.ipynb
.. _Gaia_Photometry.ipynb: Gaia_Photometry.ipynb
.. _Star_Catalog_lv0_6.ipynb: Star_Catalog_lv0_6.ipynb
.. _Star_catalog_lv7.ipynb: Star_catalog_lv7.ipynb
.. _Star_catalog_lv8.ipynb: Star_catalog_lv8.ipynb

Acknowledgement
----------------

**Infrastructure:**

* | This works has made use of the SIMBAD database, operated at CDS, Strasbourg, France.

**Software**

* | This work made use of `Astropy`_, a community-developed core Python package and an ecosystem of tools and resources for astronomy.

**Data**

* | This work made use of data from the ESA Gaia mission.
* | This work made use of data from the Hipparcos mission.


Author
-------------
-  | **Henry Leung** - henrysky_
   | University of Toronto
   | Contact Henry: henryskyleung [at] gmail.com

License
-------

The work done here is intended to be use in free and open source software (Stellarium and other open source softwares are welcome to adapt this code to generate star catalogs).

`GNU General Public License v3.0 <LICENSE>`_

.. _henrysky: https://github.com/henrysky
.. _Astropy: https://www.astropy.org
